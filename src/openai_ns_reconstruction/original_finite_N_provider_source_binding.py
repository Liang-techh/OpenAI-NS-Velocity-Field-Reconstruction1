"""Bind the original finite-``N`` C5 provider to accepted source functions.

This consumer is deliberately a source/data witness.  It authenticates the
saved reports, source recipes, and raw C5 graph archives supplied by the
checkpoint runner.  It never imports a producer or evaluates a graph node.
The generic ``C2=N*C5`` result is consumed from the saved composition report;
this module only proves that the actual first-ten C5 provider exposes the same
source handles and substitutions required by that saved algebra result.
"""

from __future__ import annotations

import ast
import base64
import gzip
import hashlib
import json
from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "original_finite_N_provider_source_binding_v1"
EXPECTED_NODE_COUNT = 205300
EXPECTED_NODES_DIGEST = "a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3"

PREFIX = "lei_ren_part1_paper_compliant_current_original_"
GENERIC_REPORT = PREFIX + "finite_N_field_density_source_composition.json"
GENERIC_REPORT_SHA256 = "0e5d3002ccdc854a1f6be32f2ea076f5c7d650ae25e6a766c7161bca4cd90e2f"
PARENT_CALCULUS_REPORT = PREFIX + "Rm_physical_Qy_Z_calculus.json"
PARENT_CALCULUS_SHA256 = "8903b59267bae01d5182bf2dcc84ed3c02c44eb8c7d3dd9d7b14af67ac092746"
MICRO_REPORT = PREFIX + "micro_evolution_source_projection.json"
MICRO_REPORT_SHA256 = "9916b2aae49a624b5b505ba0c4f7b3cfbdf92f9a7a30d7852fb247efbd085dbd"
MACRO_REPORT = PREFIX + "macro_evolution_source_projection.json"
MACRO_REPORT_SHA256 = "880abe2721144f86b0a9775f62735dccc8c5650a5fd8faab27613997daa4ae16"
R100_REPORT = PREFIX + "R100_switch_source_projection.json"
R100_REPORT_SHA256 = "92581da55a9da2a7dec6dbca76a433ce988f15786092b0ba35897b39bdcfb9d4"
R110_REPORT = PREFIX + "R110_Rm_source_projection.json"
R110_REPORT_SHA256 = "08ef44c6971d07a4c6583d481ef3a0ce1b37bf069995a732bdf2a9a7ccfa074e"

C2_PHASE = PREFIX + "C2_phase_primitives.py"
C2_PHASE_SHA256 = "84f5da3eb64096f82df88dd828b22b09bce18192a03925d91cfee0b193902d2c"
C2_PHASE_CHECK = PREFIX + "C2_phase_primitives_check.json"
C2_PHASE_CHECK_SHA256 = "76ed1af50ec799c49081ab67881b189b42e0d3137bcca770a67766d4111dff8b"
C2_PHASE_ARCHIVE = PREFIX + "C2_phase_primitives.json.gz"
C2_PHASE_ARCHIVE_SHA256 = "80e8f6df15edd94be6cec6aaf24528ab1eb1f06c7f7c5b61da47479e733a6bb0"
C2_DENSITY = PREFIX + "C2_density_transport.py"
C2_DENSITY_SHA256 = "8bfa32c87c73c5bedf3056c7171707ae934dafd950bd0df874dbd3b3450379d1"
C2_DENSITY_CHECK = PREFIX + "C2_density_transport_check.json"
C2_DENSITY_CHECK_SHA256 = "c4d2c62dd67e9dd7e171fabbdb9ba18308e420c5e0c98d01147c708840b8f583"
C2_DENSITY_ARCHIVE = PREFIX + "C2_density_transport.json.gz"
C2_DENSITY_ARCHIVE_SHA256 = "0f347824bea798ea692eca7967d7acadba49fe4336989bea1a5c86c922260877"
GENERIC_LOOP = "lei_ren_part1_paper_compliant_current_generic_loop_function_sources.py"
GENERIC_LOOP_SHA256 = "d694b8ee8a4bf1c10090d45a7f9d4d6de0792d02a1ebdd6b21872b8514a5b26f"
GENERIC_LOOP_CHECK = "lei_ren_part1_paper_compliant_current_generic_loop_function_sources_check.json"
GENERIC_LOOP_CHECK_SHA256 = "ac18d2621c67002769f8e0c76613ae91de54239699edd3801cccfaa052c19d07"
C5_SOURCE = PREFIX + "prescribed_loop_C5_functions.py"
C5_SOURCE_SHA256 = "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d"
C5_ARCHIVE = PREFIX + "prescribed_loop_C5_functions.json.gz"
C5_ARCHIVE_SHA256 = "8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa"
C5_CHECK = PREFIX + "prescribed_loop_C5_functions_check.json"
C5_CHECK_SHA256 = "2c77e80df04fae615f8ea909fc7b606bfccfa27659cb67744ae8d84020b1a9e7"
BRIDGE_SOURCE = PREFIX + "prescribed_bridge_C6_functions.py"
BRIDGE_SOURCE_SHA256 = "8ab2c26092db06453ed5a6047d8905e37fdad17eaf9c4947048498551425c605"
BRIDGE_ARCHIVE = PREFIX + "prescribed_bridge_C6_functions.json.gz"
BRIDGE_ARCHIVE_SHA256 = "5df06c4312efc01073516023fd886c071f439c11de30d3411b2e358bf76abd92"
BRIDGE_CHECK = PREFIX + "prescribed_bridge_C6_functions_check.json"
BRIDGE_CHECK_SHA256 = "9a0969a8cd4a5989dbd6272d6dd223533400113824f9352b4e57a83196bff61e"
R100_SOURCE = PREFIX + "prescribed_R100_Rm_C6_functions.py"
R100_SOURCE_SHA256 = "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32"
R100_ARCHIVE = PREFIX + "prescribed_R100_Rm_C6_functions.json.gz"
R100_ARCHIVE_SHA256 = "5e5b565dd8d2abc41e6f3df0d9e17000e524d493e58a5778dfcf217c6e42f44a"
R100_CHECK = PREFIX + "prescribed_R100_Rm_C6_functions_check.json"
R100_CHECK_SHA256 = "494caf482bffe5475c011a35a1b7d95fd02725b32eb55653f9ca25984dd98ad6"
PATCH_SOURCE = PREFIX + "true_Rm_patch_finite_N_C5_functions.py"
PATCH_SOURCE_SHA256 = "c4445ef8acfb33d29cf9022862f8b6396e6aa2f38bda7467409b5378a5deaf15"
PATCH_ARCHIVE = PREFIX + "true_Rm_patch_finite_N_C5_functions.json.gz"
PATCH_ARCHIVE_SHA256 = "8f3bc343e9dbffbdda465ce94d6183735184645a70e291cb900d3cde726769ad"
PATCH_CHECK = PREFIX + "true_Rm_patch_finite_N_C5_functions_check.json"
PATCH_CHECK_SHA256 = "efd3b9a6fdfd4522ead1806d4ef4250864ef670909ca047d0fa85d9ce23bf207"
PATCH_C6_SOURCE = PREFIX + "true_Rm_patch_C6_functions.py"
PATCH_C6_SOURCE_SHA256 = "a6f1cb5deaf17e9482aa964aef7867aadbc815e4d703a415e4945cf665f14002"
PATCH_C6_CHECK = PREFIX + "true_Rm_patch_C6_functions_check.json"
PATCH_C6_CHECK_SHA256 = "0190fcf738fb8add5be85cf0767464a39cc0d18240e6b36ad4ec7930dce38aa7"
CONTROL_CONNECTION_SOURCE = PREFIX + "whole_Z_all_N_control_connection.py"
CONTROL_CONNECTION_SOURCE_SHA256 = "aed09a89972ef9de0234e72691c4b5d33dfaa06d82d07b8df0781ef0f7a3807d"
RM_GENERIC_INPUTS_SOURCE = PREFIX + "Rm_generic_inputs.py"
RM_GENERIC_INPUTS_SOURCE_SHA256 = "62df6a6513b6e9c521cb6cf7cff25710dec17bf04f7b23f4bc0f15dbc1a57a88"
CORE_INLET_SOURCE = PREFIX + "core_inlet_axial6_generic_source.py"
CORE_INLET_SOURCE_SHA256 = "29738d46c6eb5897f0fbeefaec741d0eda4b97baed57b43ce9a7890a7dc03890"


# The runner requires every fixed record to be inherited from an accepted
# checkpoint.  Keep this table explicit: adding an unpinned path must fail at
# checkpoint preparation rather than silently broadening the source theorem.
RECORD_HASHES = {
    'lei_ren_part1_paper_compliant_current_original_switch_finite_N.py': '15a7844a1171518cd73912301407b7b7a137ecfb2f63219605283870c26e3b8d',
    'lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_correlated_functions.py': 'bb29100e671b6a1927f9897385056645a17d859108e537eec760a3904d66f832',
    'lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_source_functions.py': 'f0f6058bf77e52d9f5f04db3259ca9f737dc0fdb9281e2650b7a6ae5ca3f25fe',
    'lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_switch_functions.py': '8a3e9581a393fcbc59fdad7f24055b0176eb7edca8a76676a35c1b23b114922c',
    'lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_long_patch_functions.py': '847120150d0c8a83c2861f440c8f69514358b55633cfe8964f6de25bd5aa3ddc',

    GENERIC_REPORT: GENERIC_REPORT_SHA256,
    PARENT_CALCULUS_REPORT: PARENT_CALCULUS_SHA256,
    MICRO_REPORT: MICRO_REPORT_SHA256,
    MACRO_REPORT: MACRO_REPORT_SHA256,
    R100_REPORT: R100_REPORT_SHA256,
    R110_REPORT: R110_REPORT_SHA256,
    C2_PHASE: C2_PHASE_SHA256,
    C2_PHASE_CHECK: C2_PHASE_CHECK_SHA256,
    C2_PHASE_ARCHIVE: C2_PHASE_ARCHIVE_SHA256,
    C2_DENSITY: C2_DENSITY_SHA256,
    C2_DENSITY_CHECK: C2_DENSITY_CHECK_SHA256,
    C2_DENSITY_ARCHIVE: C2_DENSITY_ARCHIVE_SHA256,
    GENERIC_LOOP: GENERIC_LOOP_SHA256,
    GENERIC_LOOP_CHECK: GENERIC_LOOP_CHECK_SHA256,
    C5_SOURCE: C5_SOURCE_SHA256,
    C5_ARCHIVE: C5_ARCHIVE_SHA256,
    C5_CHECK: C5_CHECK_SHA256,
    BRIDGE_SOURCE: BRIDGE_SOURCE_SHA256,
    BRIDGE_ARCHIVE: BRIDGE_ARCHIVE_SHA256,
    BRIDGE_CHECK: BRIDGE_CHECK_SHA256,
    R100_SOURCE: R100_SOURCE_SHA256,
    R100_ARCHIVE: R100_ARCHIVE_SHA256,
    R100_CHECK: R100_CHECK_SHA256,
    PATCH_SOURCE: PATCH_SOURCE_SHA256,
    PATCH_ARCHIVE: PATCH_ARCHIVE_SHA256,
    PATCH_CHECK: PATCH_CHECK_SHA256,
    PATCH_C6_SOURCE: PATCH_C6_SOURCE_SHA256,
    PATCH_C6_CHECK: PATCH_C6_CHECK_SHA256,
    CONTROL_CONNECTION_SOURCE: CONTROL_CONNECTION_SOURCE_SHA256,
    RM_GENERIC_INPUTS_SOURCE: RM_GENERIC_INPUTS_SOURCE_SHA256,
    CORE_INLET_SOURCE: CORE_INLET_SOURCE_SHA256,
}

FIRST_TEN_CHARTS = (
    "first_micro", "second_micro", "frozen_macro", "first_switch", "second_switch",
    "post_power", "long_reshape", "reference", "restoration", "postrestore",
)
OUTER_SIX_CHARTS = ("Rh_reference", "O2_slope", "O2_axial", "O2_buffer", "O3_transition", "O3_power")
SIGNED_DENSITIES = ("m", "h", "k", "e", "p")
ORDINARY_Z_ORDERS = (0, 1, 2)
SOURCE_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}


def _sha(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _digest(value: Any) -> str:
    return _sha(json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False).encode("utf-8"))


def canonical_digest(value: Any) -> str:
    """Canonical digest used by the checkpoint runner."""

    return _digest(value)


def _gate(gates: list[str], name: str) -> None:
    if name not in gates:
        gates.append(name)


def _record(source_records: Mapping[str, Any], name: str, gates: list[str]) -> str | None:
    record = source_records.get(name)
    if not isinstance(record, Mapping):
        _gate(gates, f"source_record_missing:{name}")
        return None
    expected_path = "experiments/root_st073/" + name
    if record.get("path") != expected_path:
        _gate(gates, f"source_path_mismatch:{name}")
        return None
    text = record.get("source_text")
    if not isinstance(text, str):
        _gate(gates, f"source_text_missing:{name}")
        return None
    if record.get("sha256") != RECORD_HASHES[name] or _sha(text.encode("utf-8")) != RECORD_HASHES[name]:
        _gate(gates, f"source_hash_mismatch:{name}")
        return None
    return text


def _raw_archive(source_records: Mapping[str, Any], name: str, gates: list[str]) -> tuple[str | None, bytes | None]:
    """Authenticate raw bytes before accepting any decoded archive text."""

    record = source_records.get(name)
    if not isinstance(record, Mapping):
        _gate(gates, f"source_record_missing:{name}")
        return None, None
    text = record.get("source_text")
    if (record.get("path") != "experiments/root_st073/"+name
            or record.get("sha256") != RECORD_HASHES[name] or not isinstance(text, str)):
        _gate(gates, f"archive_record_contract_mismatch:{name}")
        return None, None
    encoded = record.get("raw_bytes_base64")
    if not isinstance(encoded, str):
        _gate(gates, f"raw_archive_bytes_missing:{name}")
        return None, None
    try:
        raw = base64.b64decode(encoded, validate=True)
    except (ValueError, TypeError):
        _gate(gates, f"raw_archive_base64_invalid:{name}")
        return None, None
    if _sha(raw) != RECORD_HASHES[name] or not raw.startswith(b"\x1f\x8b"):
        _gate(gates, f"raw_archive_hash_mismatch:{name}")
        return None, None
    try:
        decoded = gzip.decompress(raw)
    except (OSError, EOFError):
        _gate(gates, f"raw_archive_decode_error:{name}")
        return None, None
    if _sha(decoded) != record.get("decoded_sha256") or decoded.decode("utf-8") != text:
        _gate(gates, f"decoded_archive_mismatch:{name}")
        return None, None
    return text, raw


def _json(text: str | None, gates: list[str], name: str) -> Mapping[str, Any]:
    if text is None:
        return {}
    try:
        value = json.loads(text)
    except (TypeError, ValueError, UnicodeError):
        _gate(gates, f"json_invalid:{name}")
        return {}
    if not isinstance(value, Mapping):
        _gate(gates, f"json_not_mapping:{name}")
        return {}
    return value


def _archive(text: str | None, gates: list[str], name: str) -> Mapping[str, Any]:
    value = _json(text, gates, name)
    if not value:
        return {}
    if not isinstance(value.get("exact_graph_nodes"), list):
        _gate(gates, f"archive_graph_missing:{name}")
    return value


def _tree(text: str | None, gates: list[str], name: str) -> ast.AST:
    if text is None:
        return ast.Module(body=[], type_ignores=[])
    try:
        return ast.parse(text, filename=name)
    except (SyntaxError, ValueError):
        _gate(gates, f"source_ast_invalid:{name}")
        return ast.Module(body=[], type_ignores=[])


def _fn(tree: ast.AST, name: str) -> ast.AST | None:
    return next((node for node in ast.walk(tree) if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name), None)


def _norm(value: str) -> str:
    return "".join(value.split())


def _unparse(value: ast.AST | None) -> str:
    if value is None:
        return ""
    try:
        return ast.unparse(value)
    except (AttributeError, TypeError):
        return ""


def _has(node: ast.AST | None, *needles: str) -> bool:
    body = _norm(_unparse(node))
    return bool(node is not None and all(_norm(needle) in body for needle in needles))


def _assign_values(node: ast.AST | None, target: str) -> list[ast.AST]:
    # Parse a target as a target, so tuple parentheses are immaterial.
    try:
        expected = ast.parse(target + " = None").body[0].targets[0]
        wanted = ast.dump(expected, include_attributes=False)
    except (SyntaxError, ValueError, AttributeError):
        return []
    values: list[ast.AST] = []
    if node is None:
        return values
    for item in ast.walk(node):
        if isinstance(item, ast.Assign) and any(ast.dump(t, include_attributes=False) == wanted for t in item.targets):
            values.append(item.value)
        elif isinstance(item, ast.AnnAssign) and ast.dump(item.target, include_attributes=False) == wanted:
            values.append(item.value)
    return values


def _same_ast(value: ast.AST | None, expression: str) -> bool:
    if value is None:
        return False
    try:
        expected = ast.parse(expression, mode="eval").body
    except (SyntaxError, ValueError):
        return False
    return ast.dump(value, include_attributes=False) == ast.dump(expected, include_attributes=False)


def _assign(node: ast.AST | None, target: str, expression: str) -> bool:
    return any(_same_ast(value, expression) for value in _assign_values(node, target))


def _dict_field(node: ast.AST | None, variable: str, field: str) -> ast.AST | None:
    for value in _assign_values(node, variable):
        if isinstance(value, ast.Dict):
            for key, item in zip(value.keys, value.values):
                if isinstance(key, ast.Constant) and key.value == field:
                    return item
        if isinstance(value, ast.Call) and isinstance(value.func, ast.Name) and value.func.id == "dict":
            for keyword in value.keywords:
                if keyword.arg == field:
                    return keyword.value
    return None


def _return_field(node: ast.AST | None, field: str) -> ast.AST | None:
    """Read a concrete keyword from a function's returned ``dict`` AST."""

    if node is None:
        return None
    for item in ast.walk(node):
        if not isinstance(item, ast.Return):
            continue
        value = item.value
        if isinstance(value, ast.Call) and isinstance(value.func, ast.Name) and value.func.id == "dict":
            for keyword in value.keywords:
                if keyword.arg == field:
                    return keyword.value
        if isinstance(value, ast.Dict):
            for key, result in zip(value.keys, value.values):
                if isinstance(key, ast.Constant) and key.value == field:
                    return result
    return None


def _return_assign(node: ast.AST | None, field: str, expression: str) -> bool:
    value = _return_field(node, field)
    seen: set[str] = set()
    while isinstance(value, ast.Name) and value.id not in seen:
        seen.add(value.id)
        assignments = _assign_values(node, value.id)
        if len(assignments) != 1:
            break
        value = assignments[0]
    return _same_ast(value, expression)


def _bool(value: Any) -> bool:
    return value is True


def _projection_report(texts: Mapping[str, str], name: str, gates: list[str]) -> Mapping[str, Any]:
    report = _json(texts.get(name), gates, name)
    component = report.get("component_result") if isinstance(report.get("component_result"), Mapping) else {}
    return {"report": report, "component": component}


def _saved_background_witness(texts: Mapping[str, str], gates: list[str]) -> dict[str, Any]:
    stages = (
        (MICRO_REPORT, ("full_micro_background_evolution_projection_proved", "full_micro_own_history_projection_proved")),
        (MACRO_REPORT, ("full_macro_background_evolution_projection_proved", "full_macro_own_history_projection_proved")),
        (R100_REPORT, ("full_R100_switch_background_projection_proved", "full_R100_switch_own_history_projection_proved", "same_exact_width_source_proved")),
        (R110_REPORT, ("full_R110_Rm_background_projection_proved", "full_R110_Rm_own_history_projection_proved", "same_exact_width_source_proved")),
    )
    rows: dict[str, Any] = {}
    passed = True
    for name, flags in stages:
        item = _projection_report(texts, name, gates)
        report, component = item["report"], item["component"]
        family_ok = report.get("source_family") == SOURCE_FAMILY
        common = (
            family_ok and component.get("passed") is True and component.get("nodes_digest") == EXPECTED_NODES_DIGEST
            and component.get("node_count") == EXPECTED_NODE_COUNT and component.get("ordinary_Z_orders") == [0, 1, 2]
            and report.get("source_only") is True and report.get("numeric_values_computed") is False
            and report.get("old_operators_replayed") is False
        )
        checks = {flag: component.get(flag) is True for flag in flags}
        ok = bool(common and all(checks.values()))
        if not ok:
            _gate(gates, f"accepted_background_projection_unproved:{name}")
        rows[name] = {
            "report": name,
            "report_sha256": RECORD_HASHES.get(name),
            "schema": component.get("schema"),
            "passed": ok,
            "source_only": report.get("source_only"),
            "numeric_values_computed": report.get("numeric_values_computed"),
            "old_operators_replayed": report.get("old_operators_replayed"),
            "nodes_digest": component.get("nodes_digest"),
            "node_count": component.get("node_count"),
            "ordinary_Z_orders": component.get("ordinary_Z_orders"),
            "projection_flags": checks,
            "endpoint_ports": component.get("endpoint_ports", {}),
            "parameter_ports": component.get("parameter_ports", {}),
            "proof_premises": component.get("sublemmas", {}),
            "unresolved_gates": list(report.get("unresolved_gates", []) or component.get("unresolved_gates", []) or []),
        }
        passed = passed and ok
    parent = _json(texts.get(PARENT_CALCULUS_REPORT), gates, PARENT_CALCULUS_REPORT)
    parent_adapter_ok = parent.get("source_only") is False and parent.get("accepted_ancestors_replayed") is False
    if not parent_adapter_ok:
        _gate(gates, "frozen_calculus_parent_adapter_unproved")
    return {
        "passed": passed,
        "stages": rows,
        "frozen_calculus_parent": {
            "report": PARENT_CALCULUS_REPORT,
            "report_sha256": RECORD_HASHES.get(PARENT_CALCULUS_REPORT),
            "source_only": parent.get("source_only"),
            "accepted_ancestors_replayed": parent.get("accepted_ancestors_replayed"),
            "adapter_contract": "source_only=False is retained for the frozen calculus parent; accepted background projections are newer source-only premises",
            "passed": parent_adapter_ok,
        },
        "older_Rm_generic_input_flags_ignored": True,
    }


def _generic_saved_proof(texts: Mapping[str, str], gates: list[str]) -> dict[str, Any]:
    report = _json(texts.get(GENERIC_REPORT), gates, GENERIC_REPORT)
    component = report.get("component_result") if isinstance(report.get("component_result"), Mapping) else {}
    check = component.get("generic_density_identity_check") if isinstance(component.get("generic_density_identity_check"), Mapping) else {}
    api = check.get("public_api") if isinstance(check.get("public_api"), Mapping) else {}
    checks = api.get("checks") if isinstance(api.get("checks"), Mapping) else {}
    ok = (
        report.get("source_only") is True and report.get("numeric_values_computed") is False
        and report.get("old_operators_replayed") is False and component.get("passed") is False
        and component.get("generic_density_identity_proved") is True and check.get("proved") is True
        and api.get("proved") is True and len(checks) == 13 and all(value is True for value in checks.values())
        and component.get("nodes_digest") == EXPECTED_NODES_DIGEST and component.get("ordinary_Z_orders") == [0, 1, 2]
    )
    if not ok:
        _gate(gates, "saved_generic_density_proof_unproved")
    return {
        "passed": ok,
        "report": GENERIC_REPORT,
        "report_sha256": RECORD_HASHES.get(GENERIC_REPORT),
        "generic_density_identity_proved": component.get("generic_density_identity_proved") is True,
        "generic_density_identity_reused": ok,
        "generic_density_identity_check": check,
        "public_api_checks": dict(checks),
        "public_api_check_count": len(checks),
        "ordinary_Z_orders": component.get("ordinary_Z_orders"),
        "saved_report_passed": component.get("passed"),
        "source_only": report.get("source_only"),
        "numeric_values_computed": report.get("numeric_values_computed"),
    }


def _find_parameter(nodes: Sequence[Any], name: str) -> tuple[int, Mapping[str, Any]] | None:
    for index, row in enumerate(nodes):
        if not isinstance(row, Mapping):
            continue
        if name == "N":
            if row.get("operation") == "exact_positive_integer_expression":
                return index, row
        elif row.get("operation") == "current_original_source_parameter" and row.get("name") == name:
            return index, row
    return None


def _parameter_witness(
    c2_archive: Mapping[str, Any],
    c5_archive: Mapping[str, Any],
    r100_archive: Mapping[str, Any],
    gates: list[str],
) -> dict[str, Any]:
    """Compare the shared parameter *rows*, while qualifying every port owner."""

    c2_nodes = c2_archive.get("exact_graph_nodes") if isinstance(c2_archive.get("exact_graph_nodes"), list) else []
    c5_nodes = c5_archive.get("exact_graph_nodes") if isinstance(c5_archive.get("exact_graph_nodes"), list) else []
    r100_nodes = r100_archive.get("exact_graph_nodes") if isinstance(r100_archive.get("exact_graph_nodes"), list) else []
    shared: dict[str, Any] = {}
    shared_ok = True
    for name in ("N", "eta", "d_star", "sc", "hbB", "hbS"):
        c2 = _find_parameter(c2_nodes, name)
        c5 = _find_parameter(c5_nodes, name)
        row_ok = bool(c2 and c5 and c2[1] == c5[1])
        if row_ok:
            if name == "N":
                row_ok = c2[1].get("one_common_N_for_all_charts_and_control_denominators") is True
            else:
                row_ok = (c2[1].get("Z_independent") is True
                    and c2[1].get("quantity_is_exact_original_parameter_not_range_endpoint") is True)
        if not row_ok:
            _gate(gates, f"parameter_source_row_mismatch:{name}")
        shared_ok = shared_ok and row_ok
        shared[name] = {
            "C2": {"archive": C2_PHASE_ARCHIVE, "graph_node_count": len(c2_nodes), "index": c2[0], "row_sha256": _digest(c2[1])} if c2 else None,
            "C5": {"archive": C5_ARCHIVE, "graph_node_count": len(c5_nodes), "index": c5[0], "row_sha256": _digest(c5[1])} if c5 else None,
            "exact_row_match": row_ok,
        }
    c5_binding = c5_archive.get("prescribed_loop_C5_bindings") if isinstance(c5_archive.get("prescribed_loop_C5_bindings"), Mapping) else {}
    parameter_names = {"N": "N", "eta": "eta", "dstar": "d_star", "sc": "sc"}
    binding_ok = all(
        isinstance(c5_binding.get(binding_name), int)
        and _find_parameter(c5_nodes, source_name) is not None
        and c5_binding.get(binding_name) == _find_parameter(c5_nodes, source_name)[0]
        for binding_name, source_name in parameter_names.items()
    )
    if not binding_ok:
        _gate(gates, "C5_fixed_parameter_binding_mismatch")
    r100_binding = r100_archive.get("prescribed_R100_Rm_C6_bindings") if isinstance(r100_archive.get("prescribed_R100_Rm_C6_bindings"), Mapping) else {}
    pstar = r100_binding.get("Pstar")
    p0 = r100_binding.get("P0")
    pstar_ok = isinstance(pstar, int) and 0 <= pstar < len(r100_nodes)
    p0_ok = isinstance(p0, list) and len(p0) == 8 and all(isinstance(item, int) and 0 <= item < len(r100_nodes) for item in p0)
    if not pstar_ok:
        _gate(gates, "R100_Pstar_port_invalid")
    if not p0_ok:
        _gate(gates, "R100_P0_port_invalid")
    prefix_ok = len(c2_nodes)>=23 and len(c5_nodes)>=23 and c2_nodes[:23]==c5_nodes[:23]
    if not prefix_ok:
        _gate(gates, "common_fixed_parameter_prefix_mismatch")
    passed = bool(shared_ok and binding_ok and pstar_ok and p0_ok and prefix_ok)
    return {
        "passed": passed,
        "shared_parameter_rows": shared,
        "common_parameter_prefix": dict(node_count=23,exact_prefix_match=prefix_ok,
            C2_sha256=_digest(c2_nodes[:23]),C5_sha256=_digest(c5_nodes[:23]),
            source_values_evaluated=False),
        "C5_parameter_binding": {"archive": C5_ARCHIVE, "graph_node_count": len(c5_nodes), "N": c5_binding.get("N"), "eta": c5_binding.get("eta"), "dstar": c5_binding.get("dstar"), "sc": c5_binding.get("sc")},
        "normalized_source_unit_ports": {
            "Pstar": {"archive": R100_ARCHIVE, "graph_node_count": len(r100_nodes), "index": pstar, "row_sha256": _digest(r100_nodes[pstar]) if pstar_ok else None},
            "P0": {"archive": R100_ARCHIVE, "graph_node_count": len(r100_nodes), "indices": p0, "row_sha256": [_digest(r100_nodes[item]) for item in p0] if p0_ok else None},
        },
        "bare_indices_not_equated_across_archives": True,
    }


def _source_control_witness(trees: Mapping[str, ast.AST], texts: Mapping[str, str], gates: list[str]) -> dict[str, Any]:
    """Bind the actual root recipes, phase substitution, and inverse theorem."""

    control = trees.get(CONTROL_CONNECTION_SOURCE, ast.Module(body=[], type_ignores=[]))
    c2_phase = trees.get(C2_PHASE, ast.Module(body=[], type_ignores=[]))
    c5 = trees.get(C5_SOURCE, ast.Module(body=[], type_ignores=[]))
    bridge = trees.get(BRIDGE_SOURCE, ast.Module(body=[], type_ignores=[]))
    inputs = trees.get(RM_GENERIC_INPUTS_SOURCE, ast.Module(body=[], type_ignores=[]))
    core = trees.get(CORE_INLET_SOURCE, ast.Module(body=[], type_ignores=[]))
    r100 = trees.get(R100_SOURCE, ast.Module(body=[], type_ignores=[]))
    geometry = _fn(control, "exact_geometry")
    manifest = _fn(control, "recipe_manifest")
    original_loop = _fn(control, "original_loop")
    source_packet = _fn(c2_phase, "source_packet")
    primitive = _fn(c2_phase, "primitive_C2")
    recover = _fn(inputs, "recover_inputs")
    bindings = _fn(inputs, "source_bindings")
    generic_functions = _fn(core, "generic_functions")
    implicit = _fn(c5, "implicit_inverse")
    loop = _fn(c5, "loop_functions")
    roots = _fn(c5, "generic_roots")
    chart_offsets = {
        "first_micro": _has(geometry, "first_micro=B*(x-sc/2)"),
        "second_micro": _has(geometry, "second_micro=B*(x-sc/2)"),
        "frozen_macro": _has(geometry, "frozen_macro=Y*x+B*(2*(1-x)-sc/2)"),
        "first_switch": _has(geometry, "first_switch=Y+S*x-B*sc/2"),
        "second_switch": _has(geometry, "second_switch=Y+S*(1+x)-B*sc/2"),
        "post_power": _has(geometry, "post_power=Y+sy.log(sy.Rational(11,10))*x+2*S*(1-x)-B*sc/2"),
        "long_reshape": _has(geometry, "long_reshape=L+T*x"),
        "reference": _has(geometry, "reference=L+T+(10*(C+P)-T-8)*x"),
        "restoration": _has(geometry, "restoration=L+10*(C+P)-8+x"),
        "postrestore": _has(geometry, "postrestore=L+10*(C+P)-7+x"),
    }
    quantity_paths = {
        "E": _has(manifest, "E='original_generic_source.common_velocity_E_axial5[Z_order]'"),
        "V": _has(manifest, "V='original_generic_source.common_velocity_V_axial5[Z_order]'"),
        "a": _has(manifest, "a='original_roots.a[(0,Z_order)]'"),
        "b": _has(manifest, "b='original_full_source_quotients.actual_b_axial5[Z_order]'"),
        "t0": _has(manifest, "t0='original_roots.t0[(0,Z_order)]'"),
        "p2": _has(manifest, "p2='original_roots.p2[(0,Z_order)]'"),
        "Delta": _has(manifest, "Delta='original_full_source_quotients.actual_Delta_axial5[Z_order]'"),
    }
    inverse = {
        "C5_implicit_lambda": _has(implicit, "g.quotient(known,lam,'original phase inverse lambda=K*(1+t^2)>0')"),
        "C5_phase_endpoints": _has(loop, "angle_lower=g.zero.node", "angle_upper=twopi.node", "phase_endpoints=[0,1]"),
        "C5_unique_monotone_inverse": _has(loop, "branch_inverse_is_original_not_support_cap=True", "true_prescribed_source_C5_implicit_inverse=True"),
        "C2_same_original_inverse": _has(primitive, "inverse=source.FunctionRef(g,old['inverse'])", "same_original_inverse_not_a_new_branch=True"),
        "C2_inverse_jacobian": _has(primitive, "same strictly positive inverse Jacobian"),
        "control_inverse": _has(original_loop, "derivative_positive_certificate='K*(1+t^2)>0; a>0 and original q-flat extension'", "phase_endpoints=[0,1]", "branch_inverse_is_original_not_support_cap=True"),
        "Poisson_denominator": _has(c5, "1-2*r*cos(psi)+r^2=(cos(psi)-r)^2+sin(psi)^2>0"),
    }
    flat = {
        "C5_Delta_ge_eta_predicate": _has(loop, "shear=g.node('exact_real_comparison',operator='>=',left=roots['Delta'].value.node,right=eta.node)"),
        "control_Delta_ge_eta_predicate": _has(original_loop, "shear_flat=g.node('exact_real_comparison',operator='>=',left=roots['Delta'].value.node,right=eta.node)"),
        "C2_old_none_zero_A_B": _has(primitive, "if old is None", "A=C2Function(g.zero,g.zero,g.zero)", "B=C2Function(g.zero,g.zero,g.zero)", "exact_flat_from_original_power_admission=True"),
        "same_lazy_attributes": _has(primitive, "qold=g.nodes[old['q'][0]]", "attrs={key:value for key,value in qold.items() if key not in ('operation','active_body')}"),
        "first_micro_collar": _has(original_loop, "original_collar_extra_flat_condition='first_micro coordinate<=3*sc/4 on domain>=sc/2'"),
    }
    phase = {
        "control_Y_source": _has(geometry, "Y=sy.log(25)+1000+4*P"),
        "bridge_Y_source": _has(bridge, "self.Y=g.sub(g.unary('log',c(100)),g.unary('log',Ra))"),
        "C5_origin": _has(c5, "self.phase_origin=g.mul(background.bridge.hb,self.sc,g.constant('1/2'))"),
        "C5_offset_and_fractional_part": _has(c5, "offset=g.sub(g.unary('log',g.quotient(chart['radius'],bridge.Ra,'same positive actual R/Ra')),self.phase_origin)", "phase=g.unary('fractional_part',g.mul(self.N,offset))"),
        "C2_free_phase_packet": _has(source_packet, "source_radius_and_phase_Z_ZZ_independent=True", "same_current_Z_Taylor_coordinate=True"),
        "explicit_substitution": "phase_C2 := fractional_part(N*(log(R/Ra)-hbB*sc/2))",
    }
    normalized = {
        "generic_inputs_E_V": _has(recover, "E=raw['velocity']['theta'][0]", "V=scale(raw['velocity']['axial'][0],invS)", "C=add(E,scale(Ey,-2))", "B=scale(Vy,2)"),
        "generic_input_quotient_recipes": _has(recover, "quotient_recipes=dict(a='C/E',b='B/E',t0='-B/C',p1='inertial_theta/E',p2='inertial_axial/E'"),
        "generic_source_bindings": _has(bindings, "no_signed_sector_or_pressure_term_removed=True", "canonical_five_bases_and_directed_offsets_preserved=True"),
        "generic_p2_formula": _has(generic_functions, "p2=a.scale(a.div(a.add(axial_linear,a.scale(axial_quadratic,Pstar)),E,'same actual positive core inlet swirl E'),radius)"),
        "C5_fixed_source_units": _has(c5, "self.Pstar=background.Pstar", "self.P0=background.P0", "self.N=self.unique('exact_positive_integer_expression',None)", "self.eta=self.unique('current_original_source_parameter','eta')", "self.dstar=self.unique('current_original_source_parameter','d_star')", "self.sc=self.unique('current_original_source_parameter','sc')"),
        "R100_fixed_source_units": _has(r100, "self.F0=F0=a.trim(bridge_owner.core7['F0'])", "self.P0=bridge_owner.P0", "self.Pstar=Pstar=bridge_owner.inlet.Pstar"),
        "genuine_derivative_profiles": _has(recover, "derivative_orders=dict(profile_and_histories=5,radial_and_full_inertial=4)") and _has(core, "a=JetAlgebra(g,5)", "ordinary_slow_Z_derivative_order=k") and _has(core, "def extend6", "while len(out)<7", "A genuine order-six source is required"),
    }
    checks = {
        "all_first_ten_offsets": all(chart_offsets.values()),
        "exact_root_quantity_paths": all(quantity_paths.values()),
        "actual_C2_source_packet": _has(source_packet, "rows=dict(E=raw['common_velocity_E_axial5'],V=raw['common_velocity_V_axial5']", "projected[name]=[row[0],row[1],row[2]*2]", "range_rows_are_not_function_values=True"),
        "actual_generic_recovery": all(normalized.values()),
        "phase_substitution": all(value is True for key, value in phase.items() if key != "explicit_substitution"),
        "same_monotone_inverse": all(value is True for value in inverse.values()),
        "same_lazy_flat_branch": all(flat.values()),
        "ordinary_derivative_conversion_and_smoothness": normalized["genuine_derivative_profiles"] and _has(source_packet, "coefficient2_converted_to_ordinary_derivative_once=True") and _has(c5, "math.factorial(n)", "ordinary=FunctionJet"),
    }
    passed = bool(all(checks.values()))
    if not passed:
        _gate(gates, "actual_source_function_equality_theorem_unproved")
    return {
        "passed": passed,
        "checks": checks,
        "chart_offsets": chart_offsets,
        "root_quantity_paths": quantity_paths,
        "phase_binding": phase,
        "inverse_uniqueness_premises": inverse,
        "flat_branch_premises": flat,
        "normalized_source_units": normalized,
        "same_actual_analytic_source_domain": passed,
        "domain_contract": {
            "identity_is_conditional_on_common_original_domain": True,
            "fixed_parameters": ["N is the same positive Z-independent integer", "eta>0", "d_star>0", "Pstar>0"],
            "active_branch": ["E>0", "a>0", "Delta<eta", "gamma=2*eta-Delta>0"],
            "flat_branch": "Delta>=eta or first_micro coordinate<=3*sc/4; lazy active body is not evaluated",
            "underlying_regularities": "All-Z background function equality and genuine per-side C3 or higher source definitions imply the required root Z2 equality.",
            "numerical_admissibility_or_conditioning_proved": False,
        },
        "numerical_cone_or_admissibility_margins_claimed": False,
        "C2_free_phase_substitution": phase["explicit_substitution"],
        "C2_C5_angle_variables_alpha_renamed": {"C2": "loop_angle_<chart>", "C5": "true_loop_C5_angle_<chart>"},
    }


def _strong_phase_equivalence_witness(trees, gates):
    """Static source witness for the first-ten C2/C5 phase primitive theorem.

    This helper is intentionally source-only.  It uses the already defined
    AST helpers and checks actual operation assignments.  The two source
    graphs use different scalar-builder names, so the returned theorem steps
    state the scalar normalization used to compare them rather than claiming
    literal AST identity between C2 and C5 expressions.
    """
    empty = ast.Module(body=[], type_ignores=[])
    control = trees.get(CONTROL_CONNECTION_SOURCE, empty)
    c2 = trees.get(C2_PHASE, empty)
    c5 = trees.get(C5_SOURCE, empty)
    bridge = trees.get(BRIDGE_SOURCE, empty)
    r100 = trees.get(R100_SOURCE, empty)

    geometry = _fn(control, "exact_geometry")
    build_graph = _fn(control, "build_graph")
    c2_loop = _fn(control, "original_loop")
    c2_packet = _fn(c2, "source_packet")
    c2_primitive = _fn(c2, "primitive_C2")
    c5_loop = _fn(c5, "loop_functions")
    c5_lazy = _fn(c5, "lazy")
    c5_implicit = _fn(c5, "implicit_inverse")
    c5_roots = _fn(c5, "generic_roots")

    def require(name, value):
        ok = bool(value)
        if not ok:
            _gate(gates, "strong_phase_equivalence_unproved:" + name)
        return ok

    # The C2 and C5 scalar builders differ, but these are the actual
    # assignments that define the common normalized theorem steps.
    c2_q = require("C2_q", all((
        _assign(c2_loop, "gamma", "sub(mul(two,eta),Delta.value)"),
        _assign(c2_loop, "cut", "sub(one,div(Delta.value,eta))"),
        _assign(c2_loop, "sigma", "fun('original_flat_sigma',cut)"),
        _assign(c2_loop, "root", "fun('positive_sqrt',div(gamma,mul(two,a.value)))"),
        _assign(c2_loop, "q", "flat(mul(sigma,root))"),
    )))
    c5_q = require("C5_q", all((
        _assign(c5_loop, "gamma", "a.sub(a.fixed(g.mul(c(2),eta)),roots['Delta'])"),
        _assign(c5_loop, "cut", "a.sub(one,a.scale(roots['Delta'],g.quotient(g.one,eta,'same positive canonical eta')))"),
        _assign(c5_loop, "q", "lazy(a.mul(a.sigma(cut),a.sqrt(a.div(gamma,a.scale(roots['a'],2),'same original active a>0'))))"),
    )))
    c2_u_r_alpha = require("C2_u_r_alpha", all((
        _assign(c2_loop, "u", "div(mul(p2.value,q),dstar)"),
        _assign(c2_loop, "hinv", "div(one,fun('positive_sqrt',add(one,mul(u,u))))"),
        _assign(c2_loop, "r", "mul(u,hinv)"),
        _assign(c2_loop, "alpha", "mul(two,q,hinv)"),
    )))
    c5_u_r_alpha = require("C5_u_r_alpha", all((
        _assign(c5_loop, "u", "a.scale(a.mul(roots['p2'],q),g.quotient(g.one,dstar,'same positive canonical d_star'))"),
        _assign(c5_loop, "hinv", "a.div(one,a.sqrt(a.add(one,a.mul(u,u))),'positive sqrt(1+u^2)')"),
        _assign(c5_loop, "r", "a.mul(u,hinv)"),
        _assign(c5_loop, "alpha", "a.scale(a.mul(q,hinv),2)"),
    )))

    c2_t = require("C2_t_and_poisson", all((
        _assign(c2_loop, "nn", "sub(fun('cos',psi),r)"),
        _assign(c2_loop, "den", "add(one,neg(mul(two,r,fun('cos',psi))),mul(r,r))"),
        _assign(c2_loop, "w", "div(nn,den)"),
        _assign(c2_loop, "t", "add(t0.value,mul(alpha,w))"),
    )))
    c5_t = require("C5_t_and_poisson", all((
        _assign(c5_loop, "numerator", "m.add(cc,m.scale(rr,-1))"),
        _assign(c5_loop, "den", "m.add(m.constant(1),m.scale(m.mul(rr,cc),-2),m.mul(rr,rr))"),
        _assign(c5_loop, "t", "m.add(m.axial(roots['t0']),m.mul(m.axial(alpha),m.div(numerator,den)))"),
    )))

    c2_phi = require("C2_K_Phi", all((
        _assign(c2_loop, "ratio", "add(one,mul(t0.value,t0.value),mul(two,q,q))"),
        _assign(c2_loop, "K", "div(one,mul(twopi,ratio))"),
        _assign(c2_loop, "(T1, T2)", "(integral(t), integral(mul(t,t)))"),
        _assign(c2_loop, "P", "add(psi,T2)"),
        _assign(c2_loop, "Phi", "mul(K,P)"),
    )))
    c5_phi = require("C5_K_Phi", all((
        _assign(c5_loop, "ratio", "a.add(one,a.mul(roots['t0'],roots['t0']),a.scale(a.mul(q,q),2))"),
        _assign(c5_loop, "K", "a.div(one,a.scale(ratio,twopi),'same positive 2*pi*(1+t0^2+2*q^2)')"),
        _has(c5_loop, "T1={(k,0):integrate(t[k,0])for k in range(6)}"),
        _assign(c5_loop, "T2", "FunctionJet(tuple(g.mul(c(math.factorial(k)),integrate(t2[k,0])) for k in range(6)))"),
        _assign(c5_loop, "P", "a.add(a.fixed(angle),T2)"),
        _assign(c5_loop, "Phi", "a.mul(K,P)"),
    )))

    c2_mab = require("C2_M_A_B", all((
        _assign(c2_loop, "M", "sub(neg(div(mul(a.value,T1hat),twopi)),mul(b.value,phase))"),
        _has(c2_loop, "A=source.C1Function(flat(mul(half,a.value,chi))"),
        _has(c2_loop, "B=source.C1Function(flat(mul(half,E.value,M))"),
    )))
    c5_mab = require("C5_M_A_B", all((
        _assign(c5_loop, "M", "a.sub(a.neg(a.scale(a.mul(roots['a'],T1hat),g.quotient(g.one,twopi,'positive 2*pi'))),a.scale(roots['b'],phase))"),
        _return_assign(c5_loop, "A", "lazy(a.scale(a.mul(roots['a'],chi),'1/2'))"),
        _return_assign(c5_loop, "B_over_Pstar", "lazy(a.scale(a.mul(roots['E'],M),'1/2'))"),
    )))

    # C2 holds the phase as a source handle; the physical substitution is the
    # control graph's actual fractional-part operation, not a new symbol.
    phase_binding = require("physical_phase_substitution", all((
        _has(build_graph, "phase=g.unary('fractional_part',g.mul(N,offset))"),
        _has(c5, "phase=g.unary('fractional_part',g.mul(self.N,offset))"),
        _has(c5, "offset=g.sub(g.unary('log',g.quotient(chart['radius'],bridge.Ra,'same positive actual R/Ra')),self.phase_origin)"),
        _has(c5, "self.phase_origin=g.mul(background.bridge.hb,self.sc,g.constant('1/2'))"),
        _has(c2_packet, "source_radius_and_phase_Z_ZZ_independent=True"),
    )))

    # The ten native offsets are checked from the actual control recipe.  The
    # physical radius identities below explain why the C5 radius quotient is
    # the same offset; they also preserve hbB/hbS ownership.
    offsets = require("ten_native_offsets", all((
        _has(geometry, "first_micro=B*(x-sc/2)"),
        _has(geometry, "second_micro=B*(x-sc/2)"),
        _has(geometry, "frozen_macro=Y*x+B*(2*(1-x)-sc/2)"),
        _has(geometry, "first_switch=Y+S*x-B*sc/2"),
        _has(geometry, "second_switch=Y+S*(1+x)-B*sc/2"),
        _has(geometry, "post_power=Y+sy.log(sy.Rational(11,10))*x+2*S*(1-x)-B*sc/2"),
        _has(geometry, "long_reshape=L+T*x"),
        _has(geometry, "reference=L+T+(10*(C+P)-T-8)*x"),
        _has(geometry, "restoration=L+10*(C+P)-8+x"),
        _has(geometry, "postrestore=L+10*(C+P)-7+x"),
    )))
    radius_owner = require("hbB_hbS_radius_owners", all((
        _has(bridge, "self.Y=g.sub(g.unary('log',c(100)),g.unary('log',Ra))"),
        _has(bridge, "self.hb=hb=inlet_owner.previous.lower.parameters['hb']"),
        _has(bridge, "self.Ra=Ra=inlet_owner.radius"),
        _has(bridge, "R=g.mul(Ra,g.unary('exp',g.mul(hb,s)))"),
        _has(bridge, "Rmacro=g.mul(Ra,g.unary('exp',g.add(g.mul(c(2),hb),g.mul(length,m))))"),
        _has(r100, "self.bridge_hb=bridge_owner.hb"),
        _has(r100, "self.hb=hb=self.native_parameter('hbS'"),
        _has(r100, "R=g.mul(c(100),g.unary('exp',g.mul(hb,s)))"),
        _has(r100, "R=g.mul(c(100),g.unary('exp',g.mul(hb,g.add(g.one,t))))"),
        _has(r100, "R=g.mul(c(100),g.unary('exp',g.add(g.mul(c(2),hb),y)))"),
        _has(r100, "y=g.mul(T,q)"),
        _has(r100, "y=g.add(T,g.mul(length,q))"),
        _has(r100, "R=g.mul(c(110),g.unary('exp',y))", "R=g.mul(c(110),g.unary('exp',g.add(logref,offset)))"),
    )))

    lazy = require("shared_lazy_cutoff", all((
        _has(c2_loop, "shear_flat=g.node('exact_real_comparison',operator='>=',left=roots['Delta'].value.node,right=eta.node)"),
        _has(c5_loop, "shear=g.node('exact_real_comparison',operator='>=',left=roots['Delta'].value.node,right=eta.node)"),
        _assign(c2_loop, "collar", "g.mul(g.constant('3/4'),sc)"),
        _assign(c5_loop, "collar", "g.mul(c('3/4'),sc) if chart=='first_micro' else None"),
        _has(c2_loop, "predicate=g.node('logical_or',arguments=[shear_flat.node,collar_flat.node])"),
        _has(c5_loop, "predicate=g.node('logical_or',arguments=[shear.node,g.node('exact_real_comparison',operator='<=',left=coordinate.node,right=collar.node).node])"),
        _has(c2_loop, "active_body_not_evaluated_on_flat_branch=True"),
        _has(c5_lazy, "active_body_not_evaluated_on_flat_branch=True"),
        _has(c5_lazy, "same_original_flat_extension_through_order5=True"),
    )))

    # These are source flags/operations on both sides, not invented numeric
    # positivity certificates.  The denominator implication is mathematical:
    # r = u/sqrt(1+u^2) gives |r| < 1 for real u, hence
    # 1 - 2 r cos(psi) + r^2 = (1-r)^2 + 2 r (1-cos(psi)) > 0.
    phase_z = require("phase_Z_zero", all((
        _has(c2_loop, "phase_Z_exact_zero=True"),
        _has(c2_primitive, "phase_Z_and_ZZ_exact_zero=True"),
        _has(c2_packet, "source_radius_and_phase_Z_ZZ_independent=True"),
        _has(c5, "total_Z_phase_derivative_exact_zero=True"),
        _has(c5_loop, "angle_endpoints_independent_of_Z=True"),
    )))
    p0_owner = require("C2_exact_P0_owner", all((
        _has(c2_packet, "packet['exact_common_P0_axial5']isnotop.P0"),
        _has(c2_packet, "exact_common_P0_axial5=op.P0"),
    )))
    ordinary_rows = require("ordinary_Z_rows", all((
        _has(c2_packet, "projected[name]=[row[0],row[1],row[2]*2]"),
        _assign(c5_implicit, "ordinary", "FunctionJet((inverse,)+tuple(g.mul(g.constant(math.factorial(n)),delta[n]) for n in range(1,order+1)))"),
        _has(c2_primitive, "source_projection='same current axial coefficient[2] times 2'"),
    )))
    inverse = require("inverse_Z_and_uniqueness", all((
        _assign(c2_loop, "lam", "at(mul(K,add(one,mul(t,t))))"),
        _assign(c2_loop, "psiZ", "neg(div(at(PhiZ),lam))"),
        _has(c2_primitive, "psiZZ=g.neg(g.quotient", "g.add(at(Phi.ZZ),g.mul(g.constant(2),at(phi_psi.Z),psiZ)", "same strictly positive inverse Jacobian"),
        _assign(c5_loop, "lam", "m.mul(m.axial(K),m.add(m.constant(1),t2))"),
        _assign(c5_implicit, "lam", "at[0,1]"),
        _has(c5_implicit, "delta[n]=g.neg(g.quotient(known,lam,'original phase inverse lambda=K*(1+t^2)>0'))"),
        _has(c5, "angle_lower=g.zero.node,angle_upper=twopi.node,phase_endpoints=[0,1]"),
        _has(c2_loop, "angle_lower=g.zero.node", "angle_upper=twopi.node", "phase_endpoints=[0,1]"),
    )))
    positivity = require("common_original_admissible_domain", all((
        _has(c5_roots, "divide=lambda jet:a.div(jet,E,'same strictly positive actual prescribed swirl E')"),
        _has(c5_roots, "t0=a.neg(a.div(B,C,'same original a>0 and C=a*E'))"),
        _has(c5_roots, "a.div(a.mul(bb,bb),aa,'same original positive a')"),
        _has(c5_loop, "g.quotient(g.one,eta,'same positive canonical eta')"),
        _has(c5_loop, "g.quotient(g.one,dstar,'same positive canonical d_star')"),
        _has(c5_loop, "a.div(gamma,a.scale(roots['a'],2),'same original active a>0')"),
        _has(c5_loop, "operator='>='"),
    )))

    checks = {
        "q_u_r_alpha": c2_q and c5_q and c2_u_r_alpha and c5_u_r_alpha,
        "t_K_Phi": c2_t and c5_t and c2_phi and c5_phi,
        "M_A_B": c2_mab and c5_mab,
        "physical_phase_substitution": phase_binding and offsets and radius_owner,
        "shared_lazy_cutoff": lazy,
        "phase_Z_zero": phase_z,
        "exact_C2_P0_owner": p0_owner,
        "ordinary_Z_rows": ordinary_rows,
        "inverse_Z_and_uniqueness": inverse,
        "common_original_domain": positivity,
    }
    passed = bool(all(checks.values()))
    if not passed:
        _gate(gates, "strong_phase_equivalence_unproved")
    return {
        "passed": passed,
        "checks": checks,
        "normalized_theorem_steps": [
            "alpha-rename true_loop_C5_angle_<chart> to loop_angle_<chart>",
            "identify C2 and C5 scalar builders for q,u,r,alpha,t,K,Phi,M",
            "substitute C2 held phase by fractional_part(N*(log(R/Ra)-hbB*sc/2))",
            "use Phi_psi = K*(1+t^2) > 0 and common endpoints [0,1] for the unique inverse",
            "use ordinary source rows with C2 row[2]*2 and C5 factorial restoration",
            "conclude equality only on the common original admissible domain; no numerical cone claim",
        ],
        "loop_B_normalization": "B_loop = E*M/2; generic B = 2*V_y remains distinct",
        "phase_origin": "hbB*sc/2",
        "numeric_cone_claimed": False,
    }



def _strong_generic_equivalence_witness(trees, gates):
    """Check the two actual differential recovery recipes before identifying roots."""
    empty=ast.Module(body=[],type_ignores=[])
    recover=_fn(trees.get(RM_GENERIC_INPUTS_SOURCE,empty),'recover_inputs')
    core=trees.get(CORE_INLET_SOURCE,empty)
    generic=_fn(core,'generic_functions')
    bridge=_fn(trees.get(BRIDGE_SOURCE,empty),'generic_fields')
    phase=trees.get(C2_PHASE,empty)
    dispatch=_fn(phase,'leading_source')
    packet=_fn(phase,'source_packet')
    quotients=_fn(trees.get(PREFIX+'switch_finite_N.py',empty),'general_quotients')
    raw={
      'E':_assign(recover,'E',"raw['velocity']['theta'][0]"),
      'V':_assign(recover,'V',"scale(raw['velocity']['axial'][0],invS)"),
      'pressure':_assign(recover,'pressure','add(op.P0,p)'),
      'transport':_assign(recover,'transport','add(scale(mul(m,z),1-de),mul(dz(m),d))'),
      'Q':_assign(recover,'Q','divideL(add(scale(mul(V,z),2),scale(transport,-1)))'),
      'itl':_assign(recover,'itl','divideL(add(scale(E,-1),scale(h,1-de/2),scale(mul(dz(h),z),-(1-de)/2)))'),
      'itq':_assign(recover,'itq','divideL(add(scale(mul(k,z),2*de-1),scale(mul(dz(k),d),-1),mul(E,transport)))'),
      'izl':_assign(recover,'izl','divideL(add(scale(V,-1),scale(add(m,scale(mul(dz(m),z),-1)),(1-de)/2)))'),
      'izq':_assign(recover,'izq','divideL(add(mul(V,transport),scale(mul(e,z),2*de),scale(mul(dz(e),d),-1),scale(mul(pressure,z),2*(1+de)),scale(mul(dz(pressure),d),-1)))'),
      'C':_assign(recover,'C','add(E,scale(Ey,-2))'),
      'Broot':_assign(recover,'B','scale(Vy,2)'),
      'R':_assign(recover,'R','op.Rm_factor*x'),
      'It':_assign(recover,'It','scale(add(itl,scale(itq,S)),R)'),
      'Iz':_assign(recover,'Iz','scale(add(izl,scale(izq,S)),R)'),
      'den':_assign(recover,'den','mul(C,E)'),
      'kap':_assign(recover,'kap','add(mul(C,C),mul(B,B))'),
      'excess':_assign(recover,'excess','add(kap,scale(den,-2))'),
    }
    sectors={
      'transport':_assign(generic,'transport','a.add(a.mul(a.sub(one,de),a.mul(z,m)),a.mul(d,mz))'),
      'Q':_assign(generic,'Q','div(a.sub(a.scale(a.mul(z,V),2),transport))'),
      'theta_linear':_assign(generic,'theta_linear','div(a.sub(a.add(a.neg(E),a.mul(a.sub(one,half(de)),h)),half(a.mul(a.sub(one,de),a.mul(z,hz)))))'),
      'theta_quadratic':_assign(generic,'theta_quadratic','div(a.add(a.mul(a.sub(a.scale(de,2),one),a.mul(z,k)),a.neg(a.mul(d,kz)),a.mul(E,transport)))'),
      'axial_linear':_assign(generic,'axial_linear','div(a.add(a.neg(V),half(a.mul(a.sub(one,de),a.sub(m,a.mul(z,mz))))))'),
      'axial_quadratic':_assign(generic,'axial_quadratic','div(a.add(a.mul(V,transport),a.scale(a.mul(de,a.mul(z,e)),2),a.neg(a.mul(d,ez)),a.scale(a.mul(a.add(one,de),a.mul(z,pressure)),2),a.neg(a.mul(d,pz))))'),
      'p2':_assign(generic,'p2',"a.scale(a.div(a.add(axial_linear,a.scale(axial_quadratic,Pstar)),E,'same actual positive core inlet swirl E'),radius)"),
    }
    norm={
      'W':_assign(bridge,'W',"a.scale(F0,g.mul(g.unary('sqrt',g.mul(c(2),radius)),inverse))"),
      'FF':_assign(bridge,'FF','a.scale(a.mul(F0,F0),g.mul(radius,inverse,inverse))'),
      'm':_same_ast(_dict_field(bridge,'fields','m'),"a.scale(J['M'],inverse)"),
      'h':_same_ast(_dict_field(bridge,'fields','h'),"a.scale(a.mul(W,J['H']),'1/2')"),
      'k':_same_ast(_dict_field(bridge,'fields','k'),"a.scale(a.mul(W,J['K']),g.mul(c('1/2'),inverse))"),
      'e':_same_ast(_dict_field(bridge,'fields','e'),"a.sub(a.scale(J['A'],g.mul(inverse,inverse)),a.mul(FF,J['B']))"),
      'p':_same_ast(_dict_field(bridge,'fields','p'),"a.mul(FF,J['C'])"),
      'pressure':_assign(bridge,"fields['pressure']","a.add(fields['P0'],fields['p'])"),
      'C2_raw_normalized_m_k':_has(recover,"scale(rows[0],invS) if key in ('m','k') else rows[0]"),
      'C2_P0_same_object':_has(recover,"packet['original_P0_normalized_axial5'] is not op.P0"),
    }
    frontends={
      'early':_has(dispatch,'if chart in CHARTS[:3]', 'self.outer.owner.source_query(ends,chart,left,right)'),
      'switch':_has(dispatch,'if chart in CHARTS[3:6]','WholeZAllNSwitchFunctions.leading_packet(self.switch,ends,chart,left,right)'),
      'long':_has(dispatch,'if chart in CHARTS[6:11]','WholeZAllNLongPatchFunctions.leading_packet(self.long,ends,chart,left,right)'),
      'same_family_P0':_has(packet,"packet['source_identity']!=self.identity", "packet['exact_common_P0_axial5']isnotop.P0"),
      'queries_not_values':_has(packet,'range_rows_are_not_function_values=True','source_functions=encoded(self.functions[chart])'),
      'signed_quotients':_has(quotients,"parameters.scale(sectors[r],proxy.Pstar)","p2=parameters.scale(inertial['p2'],proxy.source_radius)","full_signed_p2_axial4=p2", "one_actual_radius_factor=proxy.source_radius"),
      'Delta_quotient':_assign(quotients,'kappa','parameters.add(f,a,parameters.quotient(f,square,a,pa))') and _assign(quotients,'Delta','parameters.add(f,kappa,[f.scalar(-2)]+[f.scalar(0)]*5)'),
    }
    checks={'raw_recovery':all(raw.values()),'C6_recovery':all(sectors.values()),'physical_normalization':all(norm.values()),'actual_frontends':all(frontends.values())}
    for key,ok in checks.items():
      if not ok:_gate(gates,'strong_generic_recovery_unproved:'+key)
    return dict(passed=all(checks.values()),checks=checks,raw_recovery=raw,C6_recovery=sectors,normalization=norm,frontends=frontends,
      common_equations={'p2':'R*(axial_linear+Pstar*axial_quadratic)/E','Delta':'(C*C+Broot*Broot)/(C*E)-2=C/E+(Broot/E)^2/(C/E)-2','pressure':'P0+p'},
      regularity_lift='The accepted all-Z background value-function identity identifies one function. Its genuine C6 definition supplies C3 regularity to that same function, so differential recovery and two ordinary Z differentiations agree. No enclosure endpoint or numerical Z3 bound is used.')


def _actual_C2_packet_bindings(archive, gates):
    packets=archive.get('actual_current_Z2_source_packets',[])
    funcs=archive.get('current_C2_phase_functions',{})
    rows=[]
    for chart in FIRST_TEN_CHARTS:
      found=[p for p in packets if isinstance(p,Mapping) and p.get('chart')==chart]
      packet=found[0] if len(found)==1 else {}
      functions=packet.get('source_functions',{})
      flags={key:packet.get(key) is True for key in ('range_rows_are_not_function_values','coefficient2_converted_to_ordinary_derivative_once','source_radius_and_phase_Z_ZZ_independent','current_C1_root_prefix_retained')}
      family=packet.get('source_identity')==packet.get('actual_source_query_identity')==SOURCE_FAMILY
      ports=isinstance(functions,Mapping) and isinstance(funcs,Mapping) and functions==funcs.get(chart)
      ok=len(found)==1 and family and ports and all(flags.values())
      if not ok:_gate(gates,'actual_C2_packet_binding_unproved:'+chart)
      rows.append(dict(chart=chart,passed=ok,source_family_matches=family,function_ports_match=ports,flags=flags,
       source_archive=C2_PHASE_ARCHIVE,packet_path='actual_current_Z2_source_packets',packet_used_for_port_ownership_only=True,
       sampled_rows_used_for_function_values=False))
    return dict(passed=len(rows)==10 and all(row['passed'] for row in rows),rows=rows)


def _explicit_phase_compositions(c2archive, c5archive, gates, *, phase_theorem_proved,
        root_theorem_proved, background_theorem_proved, parameter_theorem_proved):
    old=c2archive.get('exact_graph_nodes',[])
    functions=c2archive.get('current_C2_phase_functions',{})
    charts=c5archive.get('prescribed_loop_C5_bindings',{}).get('charts',{})
    c5count=len(c5archive.get('exact_graph_nodes',[]))
    out=[]
    for chart in FIRST_TEN_CHARTS:
      c2=functions.get(chart,{})
      c5=charts.get(chart,{})
      inverse=c2.get('inverse')
      node=old[inverse] if _valid_port(inverse,len(old),len(old)) else {}
      phase=node.get('target_phase') if isinstance(node,Mapping) else None
      premises = all(value is True for value in (phase_theorem_proved,root_theorem_proved,background_theorem_proved,parameter_theorem_proved))
      ok=(premises and node.get('operation')=='original_monotone_phase_inverse' and _valid_port(phase,len(old),len(old))
          and _valid_port(c5.get('actual_phase'),c5count,c5count)
          and _valid_port(c5.get('coordinate'),c5count,c5count) and node.get('phase_endpoints')==[0,1])
      if not ok:_gate(gates,'explicit_phase_composition_unproved:'+chart)
      out.append(dict(chart=chart,defined=ok,operation='substitute_in_source_function',
        held_phase_port=dict(archive=C2_PHASE_ARCHIVE,node=phase),
        replacement_port=dict(archive=C5_ARCHIVE,node=c5.get('actual_phase')),
        native_coordinate_port=dict(archive=C5_ARCHIVE,node=c5.get('coordinate')),
        source_A_B_jets=dict(A=c2.get('A'),B=c2.get('B')),
        root_function_correspondence={key:dict(C2=c2.get('source_roots',{}).get(key),C5=c5.get('roots',{}).get(key,[])[:3]) for key in ('E','V','a','b','t0','p2','Delta')},
        root_correspondence_basis='accepted all-Z background identity and exact generic differential recovery',
        phase_depends_only_on_native_coordinate=bool(phase_theorem_proved and parameter_theorem_proved),
        ordinary_Z_substitution_commutes=bool(premises),
        accepted_premises=dict(background_all_Z_function_identity=background_theorem_proved,
            root_recovery_and_regularity=root_theorem_proved,phase_substitution=phase_theorem_proved,
            fixed_parameter_identity=parameter_theorem_proved),
        ordinary_Z_orders=[0,1,2],naked_node_indices_are_not_cross_archive_equalities=True))
    return dict(passed=len(out)==10 and all(row['defined'] for row in out),definitions=out)

def _formula_witness(trees: Mapping[str, ast.AST], texts: Mapping[str, str], gates: list[str]) -> dict[str, Any]:
    c2 = trees.get(C2_DENSITY, ast.Module(body=[], type_ignores=[]))
    c2_phase = trees.get(C2_PHASE, ast.Module(body=[], type_ignores=[]))
    c5 = trees.get(C5_SOURCE, ast.Module(body=[], type_ignores=[]))
    bridge = trees.get(BRIDGE_SOURCE, ast.Module(body=[], type_ignores=[]))
    roots = _fn(c5, "generic_roots")
    finite = _fn(c5, "finite_N_functions")
    loop = _fn(c5, "loop_functions")
    coefficient = _fn(c2, "coefficient_C2")
    generic_fields = _fn(bridge, "generic_fields")
    source_packet = _fn(c2_phase, "source_packet")

    c2_checks = {
        "argument": _assign(coefficient, "arg", "g.mul(A.value, invN)"),
        "exprel": _assign(coefficient, "exprel", "g.unary('exprel', arg)"),
        "exponential": _assign(coefficient, "exponential", "g.unary('exp', arg)"),
        "FZZ": _assign(coefficient, "FZZ", "g.add(g.mul(E.ZZ, A.value, exprel), g.mul(g.constant(2), E.Z, A.Z, exponential), g.mul(E.value, exponential, g.add(A.ZZ, g.mul(A.Z, A.Z, invN))))"),
        "first": _same_ast(_dict_field(coefficient, "first", "m"), "B") and _same_ast(_dict_field(coefficient, "first", "h"), "F") and _same_ast(_dict_field(coefficient, "first", "k"), "alg.add(VF, EB)"),
        "second": _same_ast(_dict_field(coefficient, "second", "k"), "alg.mul(F, B)"),
        "full_second": _has(coefficient, "alg.add(B2, alg.neg(alg.scale(F2, '1/2')))", "alg.scale(F2, '1/2')"),
    }
    c5_checks = {
        "roots_E_V": _assign(roots, "(E, V)", "(a.trim(field[key]) for key in ('E', 'V'))"),
        "roots_a_b": _assign(roots, "(aa, bb)", "(divide(C), divide(B))"),
        "roots_t0": _assign(roots, "t0", "a.neg(a.div(B, C, 'same original a>0 and C=a*E'))"),
        "roots_p2": _return_assign(roots, "p2", "generic['p2']"),
        "roots_Delta": _return_assign(roots, "Delta", "a.sub(kap, a.fixed(2))"),
        "roots_B_generic": _assign(roots, "B", "a.scale(a.trim(radial['V_y']), 2)"),
        "finite_invN": _assign(finite, "invN", "g.quotient(g.one, N, 'same exact positive integer N on every source chart')"),
        "finite_A_over_N": _assign(finite, "x", "a.scale(loop['A'], invN)"),
        "finite_exp": _assign(finite, "ex", "a.exp(x)"),
        "finite_expm1": _has(finite, "FunctionJet((g.unary('expm1', x.value),) + rows(ex)[1:])"),
        "finite_E_N": _return_assign(finite, "E_N", "a.mul(E, ex)"),
        "finite_V_N": _return_assign(finite, "V_N", "a.add(V, s)"),
        "finite_density_cross_terms": _has(finite, "a.add(a.mul(V, r), a.mul(E, s), a.mul(r, s))", "a.sub(a.add(a.scale(a.mul(V, s), 2), a.mul(s, s)), a.add(a.mul(E, r), a.scale(a.mul(r, r), '1/2')))", "a.add(a.mul(E, r), a.scale(a.mul(r, r), '1/2'))"),
        "loop_phase_origin": _has(trees.get(C5_SOURCE), "self.phase_origin=g.mul(background.bridge.hb,self.sc,g.constant('1/2'))"),
        "loop_phase_offset": _has(trees.get(C5_SOURCE), "phase=g.unary('fractional_part',g.mul(self.N,offset))", "offset=g.sub(g.unary('log',g.quotient(chart['radius'],bridge.Ra,'same positive actual R/Ra')),self.phase_origin)"),
        "loop_A": _return_assign(loop, "A", "lazy(a.scale(a.mul(roots['a'], chi), '1/2'))"),
        "loop_B": _return_assign(loop, "B_over_Pstar", "lazy(a.scale(a.mul(roots['E'], M), '1/2'))"),
        "same_inverse": _has(loop, "original_monotone_phase_inverse", "derivative_positive_certificate='K*(1+t^2)>0; same original a/q-flat construction'", "branch_inverse_is_original_not_support_cap=True", "true_prescribed_source_C5_implicit_inverse=True"),
        "lazy_cutoff": _has(loop, "lazy(a.mul", "flat_predicate", "original_collar_upper"),
        "C2_phase_source_packet": _has(source_packet, "rows=dict(E=raw['common_velocity_E_axial5'],V=raw['common_velocity_V_axial5']", "projected[name]=[row[0],row[1],row[2]*2]", "coefficient2_converted_to_ordinary_derivative_once=True", "source_radius_and_phase_Z_ZZ_independent=True", "range_rows_are_not_function_values=True"),
    }
    bridge_checks = {
        "E_field": _same_ast(_dict_field(generic_fields, "fields", "E"), "a.mul(W, actual['phi'])"),
        "V_field": _same_ast(_dict_field(generic_fields, "fields", "V"), "a.scale(actual['U'], inverse)"),
        "P0_field": _same_ast(_dict_field(generic_fields, "fields", "P0"), "a.trim(P0)"),
        "E_y": _same_ast(_dict_field(generic_fields, "radial", "E_y"), "a.mul(W, a.add(a.scale(actual['phi'], '1/2'), actual['phi_y']))"),
        "V_y": _same_ast(_dict_field(generic_fields, "radial", "V_y"), "a.scale(actual['U_y'], inverse)"),
        "genericB_distinct": _assign(roots, "B", "a.scale(a.trim(radial['V_y']), 2)"),
    }
    implicit = _fn(c5, "implicit_inverse")
    positive = {
        "E_positive": _has(roots, "same strictly positive actual prescribed swirl E"),
        "a_C_positive": _has(roots, "same original a>0 and C=a*E"),
        "a_positive": _has(roots, "same original positive a"),
        "Pstar_positive": _has(generic_fields, "same strictly positive Pstar"),
        "phase_lambda_positive": _has(loop, "same positive 2*pi*(1+t0^2+2*q^2)") and _has(implicit, "original phase inverse lambda=K*(1+t^2)>0"),
        "poisson_denominator_positive": _has(c2_phase, "original Poisson denominator squared positive") and _has(loop, "den=m.add(m.constant(1),m.scale(m.mul(rr,cc),-2),m.mul(rr,rr))") and _has(c5, "1-2*r*cos(psi)+r^2=(cos(psi)-r)^2+sin(psi)^2>0 since |r|<1"),
        "IFT_unique_original_inverse": _has(loop, "true_prescribed_source_C5_implicit_inverse=True", "phase_endpoints=[0,1]", "angle_lower=g.zero.node", "angle_upper=twopi.node"),
        "flat_original_branch": _has(loop, "same original active a>0", "same positive canonical eta", "original_collar_upper"),
    }
    result = {"c2_density": c2_checks, "c5_recipe": c5_checks, "background_recovery": bridge_checks, "positive_denominator_and_IFT": positive}
    result["passed"] = bool(all(c2_checks.values()) and all(c5_checks.values()) and all(bridge_checks.values()) and all(positive.values()))
    if not result["passed"]:
        _gate(gates, "actual_source_recipe_AST_unproved")
    return result


def _receipt_witness(texts: Mapping[str, str], gates: list[str], *, include_patch: bool = False) -> dict[str, bool]:
    requirements = {
        C2_PHASE_CHECK: ("all_passed", "current_original_17_chart_source_Z2_rows_installed", "current_original_same_inverse_C2_phase_primitive_functions_installed"),
        C2_DENSITY_CHECK: ("all_passed", "current_original_17_chart_C2_signed_density_functions_installed"),
        C5_CHECK: ("all_passed", "current_original_true_prescribed_ten_chart_finite_N_density_C5_functions_installed"),
        BRIDGE_CHECK: ("all_passed", "current_original_true_prescribed_micro_macro_C6_functions_installed"),
        R100_CHECK: ("all_passed", "current_original_true_prescribed_R100_R110_C6_functions_installed"),
    }
    if include_patch:
        requirements.update({
            PATCH_CHECK: ("all_passed", "current_original_true_Rm_patch_finite_N_C5_transport_installed"),
            PATCH_C6_CHECK: ("all_passed", "current_original_true_Rm_patch_C6_functions_installed"),
        })
    result: dict[str, bool] = {}
    for name, keys in requirements.items():
        value = _json(texts.get(name), gates, name)
        ok = all(value.get(key) is True for key in keys)
        result[name] = ok
        if not ok:
            _gate(gates, f"accepted_source_receipt_unproved:{name}")
    # Generic loop check is descriptive support for the genericB port, not a
    # new algebra proof and not an admission of the old negative Rm inputs.
    generic = _json(texts.get(GENERIC_LOOP_CHECK), gates, GENERIC_LOOP_CHECK)
    result[GENERIC_LOOP_CHECK] = generic.get("all_passed") is True
    if not result[GENERIC_LOOP_CHECK]:
        _gate(gates, f"accepted_source_receipt_unproved:{GENERIC_LOOP_CHECK}")
    return result


def _prefix_witness(
    nodes: Sequence[Any],
    archive: Mapping[str, Any],
    length: int,
    name: str,
    gates: list[str],
    *,
    raw_authenticated: bool,
    gate_mismatch: bool = True,
) -> dict[str, Any]:
    saved = archive.get("exact_graph_nodes")
    ok = isinstance(saved, list) and len(saved) == length and list(nodes[:length]) == saved
    if not ok and gate_mismatch:
        _gate(gates, f"graph_prefix_mismatch:{name}")
    return {
        "archive": name,
        "archive_sha256": RECORD_HASHES.get(name),
        "saved_prefix_length": len(saved) if isinstance(saved, list) else None,
        "active_prefix_length": length,
        "exact_prefix_matches_active_graph": ok,
        "raw_archive_authenticated_before_decode": raw_authenticated,
    }


def _port(value: Any, node_count: int, *, prefix_limit: int | None = None) -> Any:
    """Retain graph function ports while rejecting malformed indices."""

    if isinstance(value, int):
        return value if 0 <= value < node_count and (prefix_limit is None or value < prefix_limit) else None
    if isinstance(value, list):
        return [_port(item, node_count, prefix_limit=prefix_limit) for item in value]
    if isinstance(value, dict):
        return {key: _port(item, node_count, prefix_limit=prefix_limit) for key, item in value.items()}
    return value


def _valid_port(value: Any, node_count: int, prefix_limit: int) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and 0 <= value < node_count and value < prefix_limit


def _jet_port(value: Any, node_count: int, prefix_limit: int, length: int = 6) -> bool:
    return isinstance(value, list) and len(value) == length and all(_valid_port(item, node_count, prefix_limit) for item in value)


def _optional_port(value: Any, node_count: int, prefix_limit: int) -> bool:
    return value is None or _valid_port(value, node_count, prefix_limit)


def _exact_chart_port_checks(row: Mapping[str, Any], node_count: int, prefix_limit: int) -> dict[str, bool]:
    roots = row.get("roots") if isinstance(row.get("roots"), Mapping) else {}
    loop = row.get("loop") if isinstance(row.get("loop"), Mapping) else {}
    finite = row.get("finite_N") if isinstance(row.get("finite_N"), Mapping) else {}
    fields = row.get("modified_fields") if isinstance(row.get("modified_fields"), Mapping) else {}
    root_names = ("E", "V", "a", "b", "p1", "p2", "t0", "Delta", "kappa")
    loop_jets = ("q", "inverse", "A", "B_over_Pstar", "A_phi", "B_over_Pstar_phi", "Phi", "phase_lambda_at_inverse", "T1_at_inverse", "direction_at_inverse", "inverse_normalized_increment")
    finite_jets = ("E_N", "V_N", "delta_E", "delta_V", "A_over_N")
    field_names = ("E", "V", "m", "h", "k", "e", "p", "P0", "pressure")
    increment_names = ("m", "h", "k", "e", "p")
    transport_names = ("m", "h", "k", "e", "p")
    partials = loop.get("phase_partial_coefficients")
    partials_ok = (
        isinstance(partials, list) and len(partials) == 21
        and all(isinstance(item, Mapping) and isinstance(item.get("Z"), int) and isinstance(item.get("angle"), int)
                and 0 <= item.get("Z") <= 5 and 0 <= item.get("angle") <= 5
                and _valid_port(item.get("node"), node_count, prefix_limit) for item in partials)
    )
    checks = {
        "roots_keys": tuple(roots) == root_names and all(_jet_port(roots.get(name), node_count, prefix_limit) for name in root_names),
        "loop_jets": all(_jet_port(loop.get(name), node_count, prefix_limit) for name in loop_jets),
        "loop_scalar_ports": all(_valid_port(loop.get(name), node_count, prefix_limit) for name in ("angle", "flat_predicate")),
        "loop_collar_port": _optional_port(loop.get("original_collar_upper"), node_count, prefix_limit),
        "loop_phase_partials": partials_ok,
        "finite_jets": all(_jet_port(finite.get(name), node_count, prefix_limit) for name in finite_jets),
        "finite_increment_jets": isinstance(finite.get("five_signed_density_increments"), Mapping) and all(
            _jet_port(finite["five_signed_density_increments"].get(name), node_count, prefix_limit) for name in increment_names
        ),
        "modified_field_jets": tuple(fields) == field_names and all(_jet_port(fields.get(name), node_count, prefix_limit) for name in field_names),
        "incoming_jets": isinstance(row.get("incoming"), Mapping) and tuple(row["incoming"]) == transport_names and all(
            _jet_port(row["incoming"].get(name), node_count, prefix_limit) for name in transport_names
        ),
        "history_jets": isinstance(row.get("history"), Mapping) and tuple(row["history"]) == transport_names and all(
            _jet_port(row["history"].get(name), node_count, prefix_limit) for name in transport_names
        ),
        "contribution_jets": isinstance(row.get("contributions"), Mapping) and tuple(row["contributions"]) == transport_names and all(
            _jet_port(row["contributions"].get(name), node_count, prefix_limit) for name in transport_names
        ),
        "memory_ports": isinstance(row.get("memory"), Mapping) and tuple(row["memory"]) == transport_names and all(
            _valid_port(row["memory"].get(name), node_count, prefix_limit) for name in transport_names
        ),
        "chart_scalar_ports": all(_valid_port(row.get(name), node_count, prefix_limit) for name in ("coordinate", "radius", "physical_measure", "global_offset", "actual_phase")),
        "chart_domain_ports": isinstance(row.get("domain"), list) and len(row["domain"]) == 2 and all(_valid_port(item, node_count, prefix_limit) for item in row["domain"]),
    }
    checks["all_required_ports"] = all(checks.values())
    return checks


def _chart_ports(
    archive: Mapping[str, Any],
    nodes: Sequence[Any],
    prefix_limit: int,
    chart: str,
    proof_premises: Mapping[str, Any],
    gates: list[str],
) -> dict[str, Any]:
    charts = archive.get("prescribed_loop_C5_bindings", {}).get("charts", {})
    row = charts.get(chart) if isinstance(charts, Mapping) else None
    if not isinstance(row, Mapping):
        _gate(gates, f"chart_port_row_missing:{chart}")
        return {"chart": chart, "present": False, "graph_function_ports": {}}
    port_checks = _exact_chart_port_checks(row, len(nodes), prefix_limit)
    if not port_checks["all_required_ports"]:
        for name, ok in port_checks.items():
            if not ok:
                _gate(gates, f"chart_port_invalid:{chart}:{name}")
    provider = {
        "archive": C5_ARCHIVE,
        "source_prefix_limit": prefix_limit,
        "roots": _port(row.get("roots", {}), len(nodes), prefix_limit=prefix_limit),
        "loop": _port(row.get("loop", {}), len(nodes), prefix_limit=prefix_limit),
        "finite_N": _port(row.get("finite_N", {}), len(nodes), prefix_limit=prefix_limit),
        "fields": _port(row.get("modified_fields", {}), len(nodes), prefix_limit=prefix_limit),
        "history": _port(row.get("history", {}), len(nodes), prefix_limit=prefix_limit),
        "incoming": _port(row.get("incoming", {}), len(nodes), prefix_limit=prefix_limit),
        "contributions": _port(row.get("contributions", {}), len(nodes), prefix_limit=prefix_limit),
        "memory": _port(row.get("memory", {}), len(nodes), prefix_limit=prefix_limit),
        "coordinate": _port(row.get("coordinate"), len(nodes), prefix_limit=prefix_limit),
        "radius": _port(row.get("radius"), len(nodes), prefix_limit=prefix_limit),
        "physical_measure": _port(row.get("physical_measure"), len(nodes), prefix_limit=prefix_limit),
        "global_offset": _port(row.get("global_offset"), len(nodes), prefix_limit=prefix_limit),
        "actual_phase": _port(row.get("actual_phase"), len(nodes), prefix_limit=prefix_limit),
    }
    return {
        "chart": chart,
        "present": port_checks["all_required_ports"],
        "source_provider": "CurrentOriginalPrescribedLoopC5Functions.source_charts -> generic_roots -> loop_functions -> finite_N_functions",
        "graph_function_ports": provider,
        "port_checks": port_checks,
        "proof_premises": dict(proof_premises),
        "endpoint_rows_are_values": False,
    }


def _c2_phase_ports(archive: Mapping[str, Any], charts: Sequence[str], gates: list[str]) -> list[dict[str, Any]]:
    nodes = archive.get("exact_graph_nodes") if isinstance(archive.get("exact_graph_nodes"), list) else []
    functions = archive.get("current_C2_phase_functions") if isinstance(archive.get("current_C2_phase_functions"), Mapping) else {}
    rows: list[dict[str, Any]] = []
    for chart in charts:
        row = functions.get(chart)
        if not isinstance(row, Mapping):
            _gate(gates, f"C2_phase_chart_missing:{chart}")
            rows.append({"chart": chart, "present": False})
            continue
        roots = row.get("source_roots") if isinstance(row.get("source_roots"), Mapping) else {}
        root_names = ("E", "V", "a", "b", "p2", "t0", "Delta")
        checks = {
            "root_jets": tuple(roots) == root_names and all(isinstance(roots.get(name), list) and len(roots[name]) == 3 and all(_valid_port(item, len(nodes), len(nodes)) for item in roots[name]) for name in root_names),
            "A_rows": isinstance(row.get("A"), list) and len(row["A"]) == 3 and all(_valid_port(item, len(nodes), len(nodes)) for item in row["A"]),
            "B_rows": isinstance(row.get("B"), list) and len(row["B"]) == 3 and all(_valid_port(item, len(nodes), len(nodes)) for item in row["B"]),
            "inverse_port": _valid_port(row.get("inverse"), len(nodes), len(nodes)),
            "psi_Z_port": _valid_port(row.get("psi_Z"), len(nodes), len(nodes)),
            "psi_ZZ_port": _valid_port(row.get("psi_ZZ"), len(nodes), len(nodes)),
            "Phi_C2_rows": isinstance(row.get("Phi_C2"), list) and len(row["Phi_C2"]) == 3 and all(_valid_port(item, len(nodes), len(nodes)) for item in row["Phi_C2"]),
            "Phi_psi_C1_rows": isinstance(row.get("Phi_psi_C1"), list) and len(row["Phi_psi_C1"]) == 3 and all(_valid_port(item, len(nodes), len(nodes)) for item in row["Phi_psi_C1"]),
            "Phi_psipsi_port": _valid_port(row.get("Phi_psipsi"), len(nodes), len(nodes)),
        }
        checks["all_required_ports"] = all(checks.values())
        if not checks["all_required_ports"]:
            for name, ok in checks.items():
                if not ok:
                    _gate(gates, f"C2_phase_port_invalid:{chart}:{name}")
        rows.append({
            "chart": chart,
            "present": checks["all_required_ports"],
            "archive": C2_PHASE_ARCHIVE,
            "archive_graph_node_count": len(nodes),
            "graph_function_ports": _port(dict(row), len(nodes), prefix_limit=len(nodes)),
            "port_checks": checks,
            "endpoint_rows_are_values": False,
            "source_flags": {key: row.get(key) for key in ("original_C1_A_B_prefix_unchanged", "original_flat_and_collar_predicates_unchanged", "same_original_inverse_not_a_new_branch", "phase_Z_and_ZZ_exact_zero") if key in row},
        })
    return rows


def _patch_witness(nodes: Sequence[Any], archive: Mapping[str, Any], source_text: str | None, gates: list[str]) -> dict[str, Any]:
    prefix = archive.get("exact_graph_nodes") if isinstance(archive.get("exact_graph_nodes"), list) else []
    prefix_ok = len(prefix) == 205147 and list(nodes[:205147]) == prefix
    source_tree = _tree(source_text, gates, PATCH_SOURCE)
    cls = next((n for n in ast.walk(source_tree) if isinstance(n, ast.ClassDef) and n.name == "CurrentOriginalTrueRmPatchFiniteNC5Functions"), None)
    checks = {
        "prefix_matches_active_graph": prefix_ok,
        "class_source_present": cls is not None,
        "source_N_eta_dstar_sc": _assign(cls, "self.N,self.eta,self.dstar,self.sc", "(self.loop.N,self.loop.eta,self.loop.dstar,self.loop.sc)"),
        "source_Pstar_P0": _assign(cls, "self.Pstar,self.P0,self.delta", "(leading.Pstar,leading.P0,leading.background.delta)"),
        "source_phase_origin": _has(cls, "self.phase_origin=self.loop.phase_origin"),
        "source_Rm_finite_functions": _has(cls, "roots=engine.generic_roots(g,bg,self.Pstar)", "loop=engine.loop_functions(g,roots,phase,self.eta,self.dstar,'Rm_patch',x,self.sc)", "finite=engine.finite_N_functions(g,roots,loop,self.N)"),
        "patch_ports_returned": _has(cls, "self.radial=radial.total_radial_functions", "self.generic=radial.generic_functions", "self.physical=radial.physical_pressure_fields"),
    }
    passed = bool(all(checks.values()))
    return {
        "source_composition_proved": False,
        "leading_patch_recipe_authenticated": passed,
        "leading_patch_identity_accepted": False,
        "patch_prefix": {"archive": PATCH_ARCHIVE, "archive_sha256": RECORD_HASHES.get(PATCH_ARCHIVE), "saved_prefix_length": len(prefix), "exact_prefix_matches_active_graph": prefix_ok},
        "graph_function_ports": _port(archive.get("Rm_patch_finite_N_C5_bindings", {}), len(nodes), prefix_limit=205147),
        "proof_premises": checks,
        "unresolved_gates": ["leading_patch_source_identity_unproved"],
        "endpoint_rows_are_values": False,
    }


def _module_hash() -> str | None:
    try:
        from pathlib import Path
        return _sha(Path(__file__).read_bytes())
    except OSError:
        return None


def derive(nodes: Sequence[Any], source_records: Mapping[str, Any], *, expected_nodes_digest: str) -> dict[str, Any]:
    """Inspect authenticated source records and bind the first-ten provider."""

    gates: list[str] = []
    patch_gates: list[str] = []
    patch_records = {PATCH_SOURCE, PATCH_ARCHIVE, PATCH_CHECK, PATCH_C6_SOURCE, PATCH_C6_CHECK}
    actual_digest = ""
    try:
        actual_digest = canonical_digest(nodes)
    except (TypeError, ValueError, OverflowError):
        _gate(gates, "nodes_digest_error")
    if expected_nodes_digest != EXPECTED_NODES_DIGEST:
        _gate(gates, "expected_nodes_digest_mismatch")
    if actual_digest != EXPECTED_NODES_DIGEST:
        _gate(gates, "nodes_digest_mismatch")
    valid_nodes = isinstance(nodes, Sequence) and not isinstance(nodes, (str, bytes, bytearray))
    if not valid_nodes or len(nodes) != EXPECTED_NODE_COUNT:
        _gate(gates, "node_count_mismatch")
    records = source_records if isinstance(source_records, Mapping) else {}
    if not isinstance(source_records, Mapping):
        _gate(gates, "source_records_not_mapping")

    texts: dict[str, str] = {}
    for name in RECORD_HASHES:
        target_gates = patch_gates if name in patch_records else gates
        if name.endswith(".gz"):
            continue
        text = _record(records, name, target_gates)
        if text is not None:
            texts[name] = text
    c2_phase_text, c2_phase_raw = _raw_archive(records, C2_PHASE_ARCHIVE, gates)
    c2_density_text, c2_density_raw = _raw_archive(records, C2_DENSITY_ARCHIVE, gates)
    c5_text, c5_raw = _raw_archive(records, C5_ARCHIVE, gates)
    patch_text, patch_raw = _raw_archive(records, PATCH_ARCHIVE, patch_gates)
    bridge_text, bridge_raw = _raw_archive(records, BRIDGE_ARCHIVE, gates)
    r100_text, r100_raw = _raw_archive(records, R100_ARCHIVE, gates)
    # Keep the authenticated decoded archive text in the same source-record
    # map used by the AST witnesses.  The runner has already authenticated the
    # raw bytes; this assignment does not decode or execute any producer.
    for name, text in (
        (C2_PHASE_ARCHIVE, c2_phase_text), (C2_DENSITY_ARCHIVE, c2_density_text),
        (C5_ARCHIVE, c5_text), (PATCH_ARCHIVE, patch_text),
        (BRIDGE_ARCHIVE, bridge_text), (R100_ARCHIVE, r100_text),
    ):
        if text is not None:
            texts[name] = text
    c2_phase_archive = _archive(c2_phase_text, gates, C2_PHASE_ARCHIVE)
    c2_density_archive = _archive(c2_density_text, gates, C2_DENSITY_ARCHIVE)
    c5_archive = _archive(c5_text, gates, C5_ARCHIVE)
    patch_archive = _archive(patch_text, patch_gates, PATCH_ARCHIVE)
    bridge_archive = _archive(bridge_text, gates, BRIDGE_ARCHIVE)
    r100_archive = _archive(r100_text, gates, R100_ARCHIVE)

    tree_names = tuple(name for name in RECORD_HASHES if name.endswith(".py"))
    trees = {name: _tree(texts.get(name), gates, name) for name in tree_names}
    generic = _generic_saved_proof(texts, gates)
    background = _saved_background_witness(texts, gates)
    formulas = _formula_witness(trees, texts, gates)
    source_control = _source_control_witness(trees, texts, gates)
    strong_phase = _strong_phase_equivalence_witness(trees, gates)
    strong_generic = _strong_generic_equivalence_witness(trees, gates)
    actual_packets = _actual_C2_packet_bindings(c2_phase_archive, gates)
    receipts = _receipt_witness(texts, gates, include_patch=False)
    patch_receipts = _receipt_witness(texts, patch_gates, include_patch=True)
    c5_prefix = (
        _prefix_witness(nodes, c5_archive, 121834, C5_ARCHIVE, gates, raw_authenticated=c5_raw is not None)
        if valid_nodes else {"exact_prefix_matches_active_graph": False, "raw_archive_authenticated_before_decode": False}
    )
    patch_prefix = (
        _prefix_witness(nodes, patch_archive, 205147, PATCH_ARCHIVE, patch_gates, raw_authenticated=patch_raw is not None, gate_mismatch=False)
        if valid_nodes else {"exact_prefix_matches_active_graph": False, "raw_archive_authenticated_before_decode": False}
    )
    parameter = _parameter_witness(c2_phase_archive, c5_archive, r100_archive, gates)
    explicit_compositions = _explicit_phase_compositions(c2_phase_archive, c5_archive, gates,
        phase_theorem_proved=strong_phase.get("passed") is True,
        root_theorem_proved=strong_generic.get("passed") is True,
        background_theorem_proved=background.get("passed") is True,
        parameter_theorem_proved=parameter.get("passed") is True)
    c2_phase_ports = _c2_phase_ports(c2_phase_archive, FIRST_TEN_CHARTS, gates)

    c5_source_family_ok = c5_archive.get("source_family") == SOURCE_FAMILY
    if not c5_source_family_ok:
        _gate(gates, "C5_archive_source_family_mismatch")
    patch_source_family_ok = patch_archive.get("source_family") == SOURCE_FAMILY
    if not patch_source_family_ok:
        _gate(patch_gates, "patch_archive_source_family_mismatch")
    c5_binding = c5_archive.get("prescribed_loop_C5_bindings") if isinstance(c5_archive.get("prescribed_loop_C5_bindings"), Mapping) else {}
    fixed_parameters = {
        key: c5_binding.get(key) for key in ("N", "eta", "dstar", "sc", "phase_origin")
    }
    fixed_parameters_ok = (
        all(isinstance(value, int) and 0 <= value < 121834 for value in fixed_parameters.values())
        and len(set(fixed_parameters.values())) == 5
        and parameter.get("passed") is True
    )
    if not fixed_parameters_ok:
        _gate(gates, "fixed_C5_parameters_unbound")

    proof_premises = {
        "same_E_V_a_b_t0_p2_Delta_source_recipes": source_control.get("passed") is True and strong_generic.get("passed") is True and background.get("passed") is True and actual_packets.get("passed") is True,
        "same_fixed_N_eta_dstar_sc_source_rows": parameter.get("passed") is True and source_control.get("normalized_source_units", {}).get("C5_fixed_source_units") is True,
        "same_Pstar_P0_normalized_source_units": source_control.get("normalized_source_units", {}).get("C5_fixed_source_units") is True and source_control.get("normalized_source_units", {}).get("R100_fixed_source_units") is True,
        "same_phase_origin_and_offset": strong_phase.get("passed") is True and parameter.get("passed") is True,
        "same_lazy_flat_collar_branch": strong_phase.get("checks", {}).get("shared_lazy_cutoff") is True,
        "same_monotone_inverse_and_IFT": strong_phase.get("checks", {}).get("inverse_Z_and_uniqueness") is True and strong_phase.get("checks", {}).get("t_K_Phi") is True,
        "C2_free_phase_explicit_substitution": source_control.get("phase_binding", {}).get("C2_free_phase_packet") is True,
        "same_A_equals_a_chi_over_2": formulas.get("c5_recipe", {}).get("loop_A") is True,
        "same_normalized_loop_B_equals_E_M_over_2": formulas.get("c5_recipe", {}).get("loop_B") is True,
        "generic_B_distinct_from_normalized_loop_B": formulas.get("background_recovery", {}).get("genericB_distinct") is True,
        "generic_density_identity_reused": generic.get("generic_density_identity_reused") is True,
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
    }
    chart_ports = [_chart_ports(c5_archive, nodes, 121834, chart, proof_premises, gates) for chart in FIRST_TEN_CHARTS] if valid_nodes else []
    charts_present = len(chart_ports) == len(FIRST_TEN_CHARTS) and all(row.get("present") for row in chart_ports)
    if not charts_present:
        _gate(gates, "first_ten_chart_ports_missing")

    patch = _patch_witness(nodes, patch_archive, texts.get(PATCH_SOURCE), patch_gates) if valid_nodes else {"source_composition_proved": False, "unresolved_gates": ["leading_patch_source_identity_unproved"]}
    c2_phase_ports_present = len(c2_phase_ports) == len(FIRST_TEN_CHARTS) and all(row.get("present") is True for row in c2_phase_ports)
    if not c2_phase_ports_present:
        _gate(gates, "first_ten_C2_phase_ports_missing")
    first_ten_ok = bool(
        actual_digest == EXPECTED_NODES_DIGEST and valid_nodes and len(nodes) == EXPECTED_NODE_COUNT
        and generic.get("passed") is True and background.get("passed") is True and formulas.get("passed") is True
        and source_control.get("passed") is True and strong_phase.get("passed") is True and all(receipts.values())
        and strong_generic.get("passed") is True and actual_packets.get("passed") is True and explicit_compositions.get("passed") is True
        and c5_prefix.get("exact_prefix_matches_active_graph") is True
        and c5_prefix.get("raw_archive_authenticated_before_decode") is True
        and charts_present and c2_phase_ports_present and fixed_parameters_ok and c5_source_family_ok
        and parameter.get("passed") is True
        and all(value is True for key, value in proof_premises.items() if key != "ordinary_Z_orders")
        and not gates
    )
    if not first_ten_ok:
        _gate(gates, "first_ten_provider_composition_unproved")
    # Patch admission remains a separate future gate even when its archive and
    # source recipe are authenticated.  The outer six never gate first-ten.
    patch_proved = False
    outer_unresolved = list(OUTER_SIX_CHARTS)
    result = {
        "schema": SCHEMA,
        "source_only": True,
        "passed": bool(first_ten_ok and not gates),
        "nodes_digest": actual_digest,
        "expected_nodes_digest": EXPECTED_NODES_DIGEST,
        "node_count": len(nodes) if valid_nodes else None,
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "ordinary_Z_not_temporal_n": True,
        "coveredcharts": list(FIRST_TEN_CHARTS) if first_ten_ok else [],
        "covered_charts": list(FIRST_TEN_CHARTS) if first_ten_ok else [],
        "first_ten_finite_N_field_source_composition_proved": first_ten_ok,
        "first_ten_density_C2_NC5_identity_proved": bool(first_ten_ok and generic.get("generic_density_identity_reused")),
        "patch_source_composition_proved": patch_proved,
        "generic_density_identity_proved": True if generic.get("generic_density_identity_proved") else False,
        "generic_density_identity_reused": generic.get("generic_density_identity_reused") is True,
        "generic_saved_proof": generic,
        "background_projection_premises": background,
        "source_recipe_theorem": formulas,
        "actual_source_equality_theorem": source_control,
        "strong_phase_equivalence_theorem": strong_phase,
        "strong_generic_recovery_theorem": strong_generic,
        "actual_C2_frontend_packet_bindings": actual_packets,
        "explicit_source_function_compositions": explicit_compositions,
        "saved_receipts": receipts,
        "patch_saved_receipts": patch_receipts,
        "fixed_source_parameters": fixed_parameters,
        "qualified_source_parameter_ports": parameter,
        "c2_phase_provider_ports": c2_phase_ports,
        "perchart_provider_ports": chart_ports,
        "patch_provider": patch,
        "archive_prefixes": {
            "tenchart_C5": c5_prefix,
            "patch_C5": patch_prefix,
            "C2_phase": {"archive": C2_PHASE_ARCHIVE, "sha256": RECORD_HASHES.get(C2_PHASE_ARCHIVE), "raw_archive_authenticated_before_decode": c2_phase_raw is not None, "graph_node_count": len(c2_phase_archive.get("exact_graph_nodes", [])) if c2_phase_archive else None},
            "C2_density": {"archive": C2_DENSITY_ARCHIVE, "sha256": RECORD_HASHES.get(C2_DENSITY_ARCHIVE), "raw_archive_authenticated_before_decode": c2_density_raw is not None, "graph_node_count": len(c2_density_archive.get("exact_graph_nodes", [])) if c2_density_archive else None},
        },
        "retained_archive_providers": {
            "bridge_C6": {"archive": BRIDGE_ARCHIVE, "sha256": RECORD_HASHES.get(BRIDGE_ARCHIVE), "prefix_nodes": len(bridge_archive.get("exact_graph_nodes", [])) if bridge_archive else None},
            "R100_Rm_C6": {"archive": R100_ARCHIVE, "sha256": RECORD_HASHES.get(R100_ARCHIVE), "prefix_nodes": len(r100_archive.get("exact_graph_nodes", [])) if r100_archive else None},
        },
        "endpoint_ports": {"endpoint_rows_are_values": False, "source_function_ports_retained": True},
        "numeric_values_computed": False,
        "old_operators_replayed": False,
        "mathematical_operator_replayed": False,
        "finite_N_history_identity": False,
        "same_density_certificate": False,
        "physical_Q_y_Z1": False,
        "global_five_moment_closure": False,
        "time_order_recursion_claimed": False,
        "timeorder_recursion_claimed": False,
        "full_goal_complete": False,
        "source_inspection_audit": {
            "source_records_only": True,
            "static_source_AST_graph_and_saved_JSON_data_inspection_only": True,
            "nodes_hashed_only": True,
            "new_source_consumer_calls": 1,
            "accepted_ancestor_operator_calls": 0,
            "old_mathematical_operator_calls": 0,
            "history_operator_calls": 0,
            "numeric_field_evaluated": False,
            "quadrature_or_norm_calls": 0,
            "recurrence_calls": 0,
            "target_calls": 0,
        },
        "accepted_source_reports": {
            "generic": {"path": GENERIC_REPORT, "sha256": RECORD_HASHES.get(GENERIC_REPORT)},
            "micro": {"path": MICRO_REPORT, "sha256": RECORD_HASHES.get(MICRO_REPORT)},
            "macro": {"path": MACRO_REPORT, "sha256": RECORD_HASHES.get(MACRO_REPORT)},
            "R100_switch": {"path": R100_REPORT, "sha256": RECORD_HASHES.get(R100_REPORT)},
            "R110_Rm": {"path": R110_REPORT, "sha256": RECORD_HASHES.get(R110_REPORT)},
        },
        "ignored_older_negative_source_flags": True,
        "unresolved_gates": sorted(set(gates)),
        "patch_unresolved_gates": sorted(set(patch_gates)) if patch_gates else list(patch.get("unresolved_gates", [])),
        "scope_limits": {
            "patch_source_composition_proved": False,
            "outer_six_unresolved": outer_unresolved,
            "history_requires_nonreset_incoming_and_physical_jacobian_once": True,
            "endpoint_rows_are_values": False,
            "generic_B_is_distinct_from_loop_B_over_Pstar": True,
            "patch_gates_are_separate_from_first_ten": True,
            "outer_six_do_not_gate_first_ten": True,
            "same_actual_analytic_source_domain_only": True,
            "numerical_cone_or_admissibility_margins_claimed": False,
        },
        "record_hashes": dict(RECORD_HASHES),
        "module_sha256": _module_hash(),
        "runtime_model_metadata": None,
    }
    return result


__all__ = ["EXPECTED_NODE_COUNT", "EXPECTED_NODES_DIGEST", "RECORD_HASHES", "canonical_digest", "derive"]
