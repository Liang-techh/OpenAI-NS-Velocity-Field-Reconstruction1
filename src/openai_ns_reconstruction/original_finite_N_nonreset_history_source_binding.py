"""Bind the first-ten finite-``N`` correction histories without resetting memory.

The accepted finite-``N`` provider report already admits the C2/C5 density
composition and the common ordinary ``Z`` orders.  This module consumes that
saved evidence and binds the actual source records to the Duhamel history
ports.  It is intentionally a source/data witness: it parses source AST and
saved JSON records, authenticates gzip bytes before decoding them, and never
imports or executes a producer, operator, checker, or numerical helper.

The first-ten result is a source-function statement.  In particular, the
``endpoint_ports`` rows are handles to functions and graph nodes, not values.
The C2 rows retain the old C1 prefixes in orders ``0, 1`` and append the
ordinary ``Z=2`` rows.  The transport audit then checks the same physical
partition, Jacobian, rate-zero memory, predecessor incoming port, and
endpoint substitution chain on every first-ten chart.
"""

from __future__ import annotations

import ast
import base64
import gzip
import hashlib
import json
from fractions import Fraction
from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "original_finite_N_nonreset_history_source_binding_v1"
EXPECTED_NODE_COUNT = 205300
EXPECTED_NODES_DIGEST = "a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3"

PREFIX = "lei_ren_part1_paper_compliant_current_original_"
PROVIDER_REPORT = PREFIX + "finite_N_provider_source_binding.json"
PROVIDER_REPORT_SHA256 = "1f23b6188abffb6d99cf00fb1274a58a2fba1100a326334a1a798d0b920f6c41"

C2_DENSITY = PREFIX + "C2_density_transport.py"
C2_DENSITY_SHA256 = "8bfa32c87c73c5bedf3056c7171707ae934dafd950bd0df874dbd3b3450379d1"
C2_DENSITY_ARCHIVE = PREFIX + "C2_density_transport.json.gz"
C2_DENSITY_ARCHIVE_SHA256 = "0f347824bea798ea692eca7967d7acadba49fe4336989bea1a5c86c922260877"
C2_DENSITY_CHECK = PREFIX + "C2_density_transport_check.json"
C2_DENSITY_CHECK_SHA256 = "c4d2c62dd67e9dd7e171fabbdb9ba18308e420c5e0c98d01147c708840b8f583"

CONTROL_CONNECTION = PREFIX + "whole_Z_all_N_control_connection.py"
CONTROL_CONNECTION_SHA256 = "aed09a89972ef9de0234e72691c4b5d33dfaa06d82d07b8df0781ef0f7a3807d"
CONTROL_ARCHIVE = PREFIX + "whole_Z_all_N_control_connection.json.gz"
CONTROL_ARCHIVE_SHA256 = "588f411e9122b9a76e6a42efc2b5754137f6af66a9cae42f8d6a7c2c59920de8"

C5_SOURCE = PREFIX + "prescribed_loop_C5_functions.py"
C5_SOURCE_SHA256 = "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d"
C5_ARCHIVE = PREFIX + "prescribed_loop_C5_functions.json.gz"
C5_ARCHIVE_SHA256 = "8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa"
C5_CHECK = PREFIX + "prescribed_loop_C5_functions_check.json"
C5_CHECK_SHA256 = "2c77e80df04fae615f8ea909fc7b606bfccfa27659cb67744ae8d84020b1a9e7"

CORE_INLET = PREFIX + "core_inlet_axial6_generic_source.py"
CORE_INLET_SHA256 = "29738d46c6eb5897f0fbeefaec741d0eda4b97baed57b43ce9a7890a7dc03890"

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


# This is deliberately explicit.  The runner rejects every record outside the
# accepted parent/source stage, so an accidental additional source cannot be
# silently admitted by this consumer.
RECORD_HASHES = {
    PROVIDER_REPORT: PROVIDER_REPORT_SHA256,
    C2_DENSITY: C2_DENSITY_SHA256,
    C2_DENSITY_ARCHIVE: C2_DENSITY_ARCHIVE_SHA256,
    C2_DENSITY_CHECK: C2_DENSITY_CHECK_SHA256,
    CONTROL_CONNECTION: CONTROL_CONNECTION_SHA256,
    CONTROL_ARCHIVE: CONTROL_ARCHIVE_SHA256,
    C5_SOURCE: C5_SOURCE_SHA256,
    C5_ARCHIVE: C5_ARCHIVE_SHA256,
    C5_CHECK: C5_CHECK_SHA256,
    CORE_INLET: CORE_INLET_SHA256,
    BRIDGE_SOURCE: BRIDGE_SOURCE_SHA256,
    BRIDGE_ARCHIVE: BRIDGE_ARCHIVE_SHA256,
    BRIDGE_CHECK: BRIDGE_CHECK_SHA256,
    R100_SOURCE: R100_SOURCE_SHA256,
    R100_ARCHIVE: R100_ARCHIVE_SHA256,
    R100_CHECK: R100_CHECK_SHA256,
    PATCH_SOURCE: PATCH_SOURCE_SHA256,
    PATCH_ARCHIVE: PATCH_ARCHIVE_SHA256,
    PATCH_CHECK: PATCH_CHECK_SHA256,
}

FIRST_TEN_CHARTS = (
    "first_micro", "second_micro", "frozen_macro", "first_switch", "second_switch",
    "post_power", "long_reshape", "reference", "restoration", "postrestore",
)
RATES = ("m", "h", "k", "e", "p")
RATE_VALUES = {"m": Fraction(1), "h": Fraction(3, 2), "k": Fraction(3, 2), "e": Fraction(1), "p": Fraction(0)}
ORDINARY_Z_ORDERS = (0, 1, 2)

SOURCE_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}

# The endpoints are source coordinates.  This manifest is also used in the
# returned report so a reviewer can see the exact native cells that were bound.
PARTITIONS = {
    "first_micro": (("sc_half",), (1, 4), (1, 2), (3, 4), (1, 1)),
    "second_micro": ((1, 1), (5, 4), (3, 2), (7, 4), (2, 1)),
    "frozen_macro": ((0, 1), (1, 2), (1, 1)),
    "first_switch": ((0, 1), (1, 2), (1, 1)),
    "second_switch": ((0, 1), (1, 2), (1, 1)),
    "post_power": ((0, 1), (1, 2), (1, 1)),
    "long_reshape": ((0, 1), (1, 2), (1, 1)),
    "reference": ((0, 1), (1, 2), (1, 1)),
    "restoration": ((0, 1), (1, 2), (1, 1)),
    "postrestore": ((0, 1), (1, 2), (1, 1)),
}

GEOMETRY = {
    "first_micro": {"domain": ("sc/2", "1"), "offset": "B*(x-sc/2)", "Jacobian": "B"},
    "second_micro": {"domain": ("1", "2"), "offset": "B*(x-sc/2)", "Jacobian": "B"},
    "frozen_macro": {"domain": ("0", "1"), "offset": "Y*x+B*(2*(1-x)-sc/2)", "Jacobian": "Y-2*B"},
    "first_switch": {"domain": ("0", "1"), "offset": "Y+S*x-B*sc/2", "Jacobian": "S"},
    "second_switch": {"domain": ("0", "1"), "offset": "Y+S*(1+x)-B*sc/2", "Jacobian": "S"},
    "post_power": {"domain": ("0", "1"), "offset": "Y+sy.log(sy.Rational(11,10))*x+2*S*(1-x)-B*sc/2", "Jacobian": "sy.log(sy.Rational(11,10))-2*S"},
    "long_reshape": {"domain": ("0", "1"), "offset": "L+T*x", "Jacobian": "T"},
    "reference": {"domain": ("0", "1"), "offset": "L+T+(10*(C+P)-T-8)*x", "Jacobian": "10*(C+P)-T-8"},
    "restoration": {"domain": ("0", "1"), "offset": "L+10*(C+P)-8+x", "Jacobian": "1"},
    "postrestore": {"domain": ("0", "1"), "offset": "L+10*(C+P)-7+x", "Jacobian": "1"},
}

POSITIVE_WIDTH_CONDITIONS = {
    "first_micro_and_second_micro": "0 < sc < 1/4 and hbB > 0",
    "frozen_macro": "Y-2*hbB > 0",
    "first_switch_and_second_switch": "hbS > 0",
    "post_power": "log(11/10)-2*hbS > 0",
    "long_reshape": "T > 0",
    "reference": "10*(logC+logP)-T-8 > 0",
    "restoration_and_postrestore": "1 > 0",
}


def _sha(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def canonical_digest(value: Any) -> str:
    """Return the runner's canonical JSON digest without evaluating graph nodes."""

    return _sha(json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8"))


def _gate(gates: list[str], name: str) -> None:
    if name not in gates:
        gates.append(name)


def _record(source_records: Mapping[str, Any], name: str, gates: list[str]) -> str | None:
    record = source_records.get(name)
    if not isinstance(record, Mapping):
        _gate(gates, f"source_record_missing:{name}")
        return None
    expected_path = "experiments/root_st073/" + name
    text = record.get("source_text")
    ok = (
        record.get("path") == expected_path
        and isinstance(text, str)
        and record.get("sha256") == RECORD_HASHES[name]
        and _sha(text.encode("utf-8")) == RECORD_HASHES[name]
    )
    if not ok:
        _gate(gates, f"source_record_contract_mismatch:{name}")
        return None
    return text


def _raw_archive(source_records: Mapping[str, Any], name: str, gates: list[str]) -> str | None:
    """Authenticate raw gzip bytes before accepting decoded JSON text."""

    record = source_records.get(name)
    if not isinstance(record, Mapping):
        _gate(gates, f"source_record_missing:{name}")
        return None
    text = record.get("source_text")
    encoded = record.get("raw_bytes_base64")
    if (
        record.get("path") != "experiments/root_st073/" + name
        or record.get("sha256") != RECORD_HASHES[name]
        or record.get("encoding") != "gzip"
        or not isinstance(text, str)
        or not isinstance(encoded, str)
    ):
        _gate(gates, f"archive_record_contract_mismatch:{name}")
        return None
    try:
        raw = base64.b64decode(encoded, validate=True)
    except (TypeError, ValueError):
        _gate(gates, f"raw_archive_base64_invalid:{name}")
        return None
    # The raw compressed bytes are pinned before decompression by design.
    if _sha(raw) != RECORD_HASHES[name] or not raw.startswith(b"\x1f\x8b"):
        _gate(gates, f"raw_archive_hash_mismatch:{name}")
        return None
    try:
        decoded = gzip.decompress(raw)
        decoded_text = decoded.decode("utf-8")
    except (OSError, EOFError, UnicodeError):
        _gate(gates, f"raw_archive_decode_error:{name}")
        return None
    if _sha(decoded) != record.get("decoded_sha256") or decoded_text != text:
        _gate(gates, f"decoded_archive_mismatch:{name}")
        return None
    return text


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


def _tree(text: str | None, name: str, gates: list[str]) -> ast.AST:
    if text is None:
        return ast.Module(body=[], type_ignores=[])
    try:
        return ast.parse(text, filename=name)
    except (SyntaxError, ValueError):
        _gate(gates, f"source_ast_invalid:{name}")
        return ast.Module(body=[], type_ignores=[])


def _function(tree: ast.AST, name: str) -> ast.AST | None:
    return next(
        (
            item for item in ast.walk(tree)
            if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)) and item.name == name
        ),
        None,
    )


def _method(tree: ast.AST, class_name: str, name: str) -> ast.AST | None:
    """Find a method on one named class, avoiding same-named helper methods."""

    for item in ast.walk(tree):
        if not isinstance(item, ast.ClassDef) or item.name != class_name:
            continue
        return next(
            (
                member for member in item.body
                if isinstance(member, (ast.FunctionDef, ast.AsyncFunctionDef))
                and member.name == name
            ),
            None,
        )
    return None


def _func_path(node: ast.AST) -> str:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        parent = _func_path(node.value)
        return f"{parent}.{node.attr}" if parent else node.attr
    return ""


def _literal(node: ast.AST | None) -> Any:
    if isinstance(node, ast.Constant):
        return node.value
    return None


def _same_ast(node: ast.AST | None, expression: str) -> bool:
    if node is None:
        return False
    try:
        expected = ast.parse(expression, mode="eval").body
    except (SyntaxError, ValueError):
        return False
    return ast.dump(node, include_attributes=False) == ast.dump(
        expected, include_attributes=False,
    )


def _subtree(node: ast.AST | None, expression: str) -> bool:
    if node is None:
        return False
    try:
        expected = ast.parse(expression, mode="eval").body
    except (SyntaxError, ValueError):
        return False
    wanted = ast.dump(expected, include_attributes=False)
    return any(ast.dump(item, include_attributes=False) == wanted for item in ast.walk(node))


def _assignment_matches(node: ast.AST | None, target: str, expression: str) -> bool:
    """Match an assignment target/value pair structurally."""

    try:
        expected_target = ast.parse(target + " = None").body[0].targets[0]
    except (SyntaxError, ValueError):
        return False
    target_dump = ast.dump(expected_target, include_attributes=False)
    for item in ast.walk(node) if node is not None else ():
        if isinstance(item, ast.Assign):
            targets = item.targets
            value = item.value
        elif isinstance(item, ast.AnnAssign):
            targets = [item.target]
            value = item.value
        else:
            continue
        if any(ast.dump(target_item, include_attributes=False) == target_dump for target_item in targets):
            if _same_ast(value, expression):
                return True
    return False


def _call_positional_prefix(node: ast.AST | None, path: str, expressions: Sequence[str]) -> bool:
    """Match the leading positional arguments of a call by AST."""

    expected: list[ast.AST] = []
    try:
        expected = [ast.parse(value, mode="eval").body for value in expressions]
    except (SyntaxError, ValueError):
        return False
    expected_dumps = [ast.dump(value, include_attributes=False) for value in expected]
    return any(
        len(call.args) >= len(expected_dumps)
        and all(
            ast.dump(call.args[index], include_attributes=False) == expected_dumps[index]
            for index in range(len(expected_dumps))
        )
        for call in _calls(node, path)
    )


def _calls(node: ast.AST | None, path: str) -> list[ast.Call]:
    if node is None:
        return []
    return [
        item for item in ast.walk(node)
        if isinstance(item, ast.Call) and _func_path(item.func) == path
    ]


def _call_with_keywords(
    node: ast.AST | None,
    path: str,
    required: Mapping[str, Any] | None = None,
) -> bool:
    required = required or {}
    for call in _calls(node, path):
        keywords = {item.arg: item.value for item in call.keywords if item.arg is not None}
        if all(
            key in keywords and (
                value is _MISSING or _literal(keywords[key]) == value
            )
            for key, value in required.items()
        ):
            return True
    return False


_MISSING = object()


def _assignment(node: ast.AST | None, target: str) -> list[ast.AST]:
    if node is None:
        return []
    result: list[ast.AST] = []
    try:
        expected = ast.dump(ast.parse(target + " = None").body[0].targets[0], include_attributes=False)
    except (SyntaxError, AttributeError, ValueError):
        return result
    for item in ast.walk(node):
        if isinstance(item, ast.Assign):
            targets = item.targets
        elif isinstance(item, ast.AnnAssign):
            targets = [item.target]
        else:
            continue
        if any(ast.dump(t, include_attributes=False) == expected for t in targets):
            result.append(item.value)
    return result


def _dict_item(node: ast.AST | None, target: str, key: str) -> ast.AST | None:
    for value in _assignment(node, target):
        if isinstance(value, ast.Dict):
            for item_key, item_value in zip(value.keys, value.values):
                if _literal(item_key) == key:
                    return item_value
        elif isinstance(value, ast.Call) and _func_path(value.func) == "dict":
            for keyword in value.keywords:
                if keyword.arg == key:
                    return keyword.value
    return None


def _dict_update_item(node: ast.AST | None, target: str, key: str) -> ast.AST | None:
    for call in _calls(node, f"{target}.update"):
        for keyword in call.keywords:
            if keyword.arg == key:
                value = keyword.value
                if isinstance(value, ast.Tuple) and len(value.elts) == 2:
                    return value
                return value
    return None


def _source_ast_witness(trees: Mapping[str, ast.AST], gates: list[str], patch_gates: list[str]) -> dict[str, Any]:
    """Inspect the defining source recipes with structural AST checks."""

    c2 = trees.get(C2_DENSITY)
    control = trees.get(CONTROL_CONNECTION)
    c5 = trees.get(C5_SOURCE)
    core = trees.get(CORE_INLET)
    bridge = trees.get(BRIDGE_SOURCE)
    r100 = trees.get(R100_SOURCE)
    patch = trees.get(PATCH_SOURCE)

    c2_build = _function(c2, "build")
    c2_coeff = _function(c2, "coefficient_C2")
    c5_transport = _function(c5, "transport")
    c5_finite = _function(c5, "finite_N_functions")
    c5_loop_init = _method(c5, "CurrentOriginalPrescribedLoopC5Functions", "__init__")
    control_geometry = _function(control, "exact_geometry")
    control_graph = _function(control, "build_graph")
    core_at = _function(core, "at")
    core_integral = _function(core, "integral")
    bridge_histories = _function(bridge, "histories")
    bridge_init = _method(bridge, "CurrentOriginalPrescribedBridgeC6Functions", "__init__")
    r100_init = _method(r100, "CurrentOriginalPrescribedR100RmC6Functions", "__init__")
    r100_add = _function(r100, "add_chart")
    patch_transport = _function(patch, "transport_log")
    patch_init = _method(patch, "CurrentOriginalTrueRmPatchFiniteNC5Functions", "__init__")

    # C2 appends the genuine ZZ row to each old C1 source window and integrates
    # the same source density over every native partition with one J.
    c2_checks = {
        "build_function_present": c2_build is not None,
        "coefficient_function_present": c2_coeff is not None,
        "C2_coefficient_uses_exprel_and_full_exp": (
            _subtree(c2_coeff, "g.unary('exprel',arg)")
            and _subtree(c2_coeff, "g.unary('exp',arg)")
        ),
        "C2_second_order_mixed_products": (
            _subtree(c2_coeff, "alg.mul(F,B)")
            and _subtree(c2_coeff, "alg.add(B2,alg.neg(alg.scale(F2,'1/2')))")
        ),
        "C2_rate_kernel_uses_right_minus_offset": _subtree(
            c2_build,
            "g.unary('exp',g.neg(g.mul(g.constant(rate),g.sub(ref(old['right_offset']),ref(old['offset'])))))",
        ),
        "C2_kernel_variable_is_that_definition": _assignment_matches(c2_build, "kernel",
            "g.unary('exp',g.neg(g.mul(g.constant(rate),g.sub(ref(old['right_offset']),ref(old['offset'])))))"),
        "C2_q_variable_is_actual_density_ZZ": _assignment_matches(c2_build, "q", "density[key].ZZ"),
        "C2_integrates_density_times_kernel_times_J": _subtree(
            c2_build,
            "g.mul(kernel,q,ref(old['Jacobian']))",
        ),
        "C2_integral_has_one_J_source_measure": _call_with_keywords(
            c2_build,
            "current.controls.integral",
            {
                "measure": "original dlogR: Jacobian exactly once",
                "N_scaled_density": True,
                "ordinary_slow_Z_derivative_order": 2,
                "integration_endpoints_kernel_radius_and_global_phase_Z_independent": True,
            },
        ),
        "C2_appends_ZZ_to_old_windows": (
            _subtree(c2_build, "phase.C2Function(*map(ref,old['density'][key]),q.ZZ)")
            and _subtree(c2_build, "phase.C2Function(*map(ref,old['contributions'][key]),zz)")
            and _subtree(c2_build, "g.add(g.mul(ref(old['memory'][key]),before[key].ZZ),zz)")
        ),
        "C2_retains_predecessor_ZZ_memory": (
            _call_with_keywords(c2_build, "dict", {
                "predecessor_ZZ_memory_retained": True,
                "quiet_local_zero_does_not_reset_incoming": True,
            })
        ),
    }

    # The saved control graph owns the source-defined offsets, domains, phase,
    # and physical Jacobians.  Exact expressions are compared as ASTs rather
    # than by whitespace or source-text substrings.
    offsets = _dict_item(control_geometry, "offsets", "first_micro")
    baseline_domains = _assignment(control_geometry, "domains")
    expected_baseline = ast.parse(
        "{key:(sy.Integer(0),sy.Integer(1)) for key in CHARTS}", mode="eval",
    ).body
    domain_checks = {
        "baseline_domains_zero_one": any(
            ast.dump(value, include_attributes=False)
            == ast.dump(expected_baseline, include_attributes=False)
            for value in baseline_domains
        ),
        "first_micro_domain_sc_half_one": _same_ast(
            _dict_update_item(control_geometry, "domains", "first_micro"),
            "(sc/2,sy.Integer(1))",
        ),
        "second_micro_domain_one_two": _same_ast(
            _dict_update_item(control_geometry, "domains", "second_micro"),
            "(sy.Integer(1),sy.Integer(2))",
        ),
    }
    offset_checks: dict[str, bool] = {}
    for chart, recipe in GEOMETRY.items():
        value = _dict_item(control_geometry, "offsets", chart)
        offset_checks[chart] = _same_ast(value, recipe["offset"])
    geometry_checks = {
        "exact_geometry_function_present": control_geometry is not None,
        "offsets": offset_checks,
        "offsets_all_match": bool(offset_checks) and all(offset_checks.values()),
        "domains": domain_checks,
        "control_graph_function_present": control_graph is not None,
        "all_N_phase_is_fractional_part_N_offset": _subtree(
            control_graph, "g.unary('fractional_part',g.mul(N,offset))",
        ),
        "control_integral_has_physical_J_once": _call_with_keywords(
            control_graph,
            "controls.integral",
            {
                "measure": "original dlogR: Jacobian exactly once",
                "N_scaled_density": True,
                "Z_derivative_order_held_phase": True,
            },
        ),
        "Jacobian_is_symbolic_offset_derivative": _subtree(
            control_graph,
            "sy.diff(offsets[chart],x)",
        ),
        "Jacobian_assignment_is_symbolic_offset_derivative": _assignment_matches(control_graph, "jac",
            "source.expression(g,sy.diff(offsets[chart],x),parameters,t)"),
        "control_history_is_nonreset_c1_add": _subtree(
            control_graph,
            "g.c1add(g.c1scale(memory,before[key]),contribution)",
        ),
        "control_window_keeps_source_cells": _call_with_keywords(
            control_graph,
            "dict",
            {"actual_current_source_cells": _MISSING},
        ),
    }

    # The common core substitution/integral API is source-bound.  These checks
    # ensure that endpoint substitutions are functions of the original jet,
    # not interval endpoints promoted to field values.
    core_checks = {
        "function_substitution_is_Z_independent": _call_with_keywords(
            core_at, "self.g.node", {"Z_independent_substitution": True},
        ),
        "definite_integral_uses_variable_name": _call_positional_prefix(
            core_integral,
            "controls.integral",
            ("g", "q", "g.nodes[variable.node]['name']", "lower", "upper"),
        ),
        "definite_integral_has_ordinary_Z_order": _call_with_keywords(
            core_integral,
            "controls.integral",
            {
                "endpoint_independent_of_Z": True,
                "true_core_inlet_complete_history": True,
                "numerical_cap_not_defining_integrand": True,
            },
        ),
    }

    c5_checks = {
        "same_five_exact_rate_values": _assignment_matches(c5, "RATES",
            "dict(m=Fraction(1),h=Fraction(3,2),k=Fraction(3,2),e=Fraction(1),p=Fraction(0))"),
        "transport_function_present": c5_transport is not None,
        "transport_kernel_has_physical_jacobian": _subtree(
            c5_transport,
            "g.mul(g.constant(-rate),jacobian,g.sub(coordinate,t))",
        ),
        "transport_body_has_jacobian_once": _subtree(
            c5_transport,
            "g.mul(jacobian,kernel)",
        ),
        "transport_integral_has_nonreset_receipt": _call_with_keywords(
            c5_transport,
            "a.integral",
            {
                "actual_finite_N_signed_correction": True,
                "physical_Jacobian_applied_once": True,
                "zero_correction_inlet_only_at_original_rminus": True,
                "predecessor_memory_never_reset": True,
            },
        ),
        "transport_has_rate_zero_memory_branch": _subtree(
            c5_transport,
            "g.one",
        ),
        "transport_outgoing_is_memory_plus_contribution": _subtree(
            c5_transport,
            "a.add(a.scale(incoming[name],decay),contribution)",
        ),
        "finite_N_retains_expm1_and_cross_terms": (
            _subtree(c5_finite, "g.unary('expm1',x.value)")
            and _subtree(c5_finite, "a.mul(r,s)")
            and _subtree(c5_finite, "a.scale(a.mul(r,r),'1/2')")
        ),
        "C5_endpoint_substitution_is_function_port": _subtree(
            c5_loop_init,
            "a.at(jet,coordinate,upper)",
        ),
        "C5_endpoint_becomes_next_incoming": _assignment_matches(
            c5_loop_init,
            "incoming",
            "self.endpoints[name]",
        ),
        "C5_postrestore_outgoing_is_Rm_correction": _assignment_matches(
            c5_loop_init,
            "self.Rm_correction",
            "incoming",
        ),
        "C5_chart_retains_physical_measure_source": _subtree(
            c5_loop_init,
            "chart['physical_measure']",
        ),
    }

    bridge_checks = {
        "bridge_history_function_present": bridge_histories is not None,
        "bridge_history_keeps_initial_memory": _subtree(
            bridge_histories,
            "a.scale(initial[name],decay)",
        ),
        "bridge_source_measure_once": _call_with_keywords(
            bridge_histories,
            "a.integral",
            {
                "original_incoming_memory_retained": True,
                "source_measure_applied_once": True,
            },
        ),
    }
    r100_checks = {
        "R100_chart_constructor_present": r100_add is not None,
        "R100_chart_keeps_incoming_moments": _call_with_keywords(
            r100_add,
            "dict",
            {"incoming_histories": _MISSING},
        ),
        "R100_chart_keeps_physical_measure": _call_with_keywords(
            r100_add,
            "dict",
            {"physical_measure": _MISSING},
        ),
    }
    patch_checks = {
        "patch_transport_function_present": patch_transport is not None,
        "patch_transport_uses_dt_over_t": _subtree(
            patch_transport,
            "g.quotient(g.one,t,'positive original patch integration coordinate')",
        ),
        "patch_transport_keeps_nonzero_inlet": _call_with_keywords(
            patch_transport,
            "a.integral",
            {
                "physical_Jacobian_applied_once": True,
                "nonzero_original_Rm_correction_inlet_retained": True,
                "predecessor_memory_never_reset": True,
            },
        ),
        "patch_keeps_C5_Rm_inlet": _subtree(
            patch_init,
            "self.loop.Rm_correction",
        ),
        "patch_does_not_reset_at_Rm": _call_with_keywords(
            patch_init,
            "dict",
            {"nonzero_true_Rm_correction_and_signed_pressure_inlets_retained": True},
        ),
    }

    measure_checks = {
        "micro_measure_hbB": _same_ast(_dict_item(bridge_init, "self.micro", "physical_measure"), "hb"),
        "macro_measure_length": _same_ast(_dict_item(bridge_init, "self.macro", "physical_measure"), "length"),
        "macro_length_Y_minus_two_hbB": _assignment_matches(bridge_init, "self.length", "g.sub(self.Y,g.mul(c(2),hb))"),
        "post_power_J_log_ratio_minus_two_hbS": _assignment_matches(r100_init, "self.post_length", "g.sub(g.unary('log',c('11/10')),g.mul(c(2),hb))"),
        "reference_J_complete_finite_length": _assignment_matches(r100_init, "self.reference_length", "g.sub(g.sub(logref,T),c(8))"),
        "reference_logref_same_original_parameters": _assignment_matches(r100_init, "self.logref", "g.mul(c(10),g.add(self.core.parameters['logC'],self.core.parameters['logP']))"),
        "R100_constructor_stores_its_J_argument": _same_ast(_dict_item(r100_add, "self.charts[name]", "physical_measure"), "jacobian"),
        "C5_stores_same_source_chart_J": _same_ast(_dict_item(c5_loop_init, "self.charts[name]", "physical_measure"), "chart['physical_measure']"),
    }
    native_chart_args = {
        "first_switch": ("'first_switch'", "s", "R", "hb"),
        "second_switch": ("'second_switch'", "t", "R", "hb"),
        "post_power": ("'post_switch'", "p", "R", "self.post_length"),
        "long_reshape": ("'long_reshape'", "q", "R", "T"),
        "reference": ("'reference'", "q", "R", "length"),
        "restoration": ("'restoration'", "t", "R", "g.one"),
        "postrestore": ("'post_restore'", "p", "R", "g.one"),
    }
    for chart, args in native_chart_args.items():
        measure_checks[chart+"_native_J_argument"] = _call_positional_prefix(r100_init, "self.add_chart", args)
    source_chart_values = {"first_micro": "bridge.micro", "second_micro": "bridge.micro", "frozen_macro": "bridge.macro"}
    source_chart_values.update({name:"background.charts["+repr(alias)+"]" for name,alias in (
        ("first_switch","first_switch"),("second_switch","second_switch"),("post_power","post_switch"),
        ("long_reshape","long_reshape"),("reference","reference"),("restoration","restoration"),("postrestore","post_restore"))})
    for chart, expression in source_chart_values.items():
        measure_checks[chart+"_C5_source_chart"] = _same_ast(_dict_item(c5_loop_init, "self.source_charts", chart), expression)

    groups = {
        "C2_density": c2_checks,
        "control_geometry_and_graph": geometry_checks,
        "core_substitution_and_integral_api": core_checks,
        "C5_transport_and_endpoints": c5_checks,
        "bridge_history": bridge_checks,
        "R100_history": r100_checks,
        "patch_history": patch_checks,
        "same_C5_native_physical_measure_sources": measure_checks,
    }
    for group_name, checks in groups.items():
        for check_name, passed in checks.items():
            if isinstance(passed, Mapping):
                passed = all(
                    item if isinstance(item, bool) else all(item.values())
                    for item in passed.values()
                )
            if not passed:
                _gate(patch_gates if group_name == "patch_history" else gates,
                    f"source_ast_unproved:{group_name}:{check_name}")
    return groups


def _check_receipt(
    report: Mapping[str, Any],
    required_flags: Sequence[str],
    name: str,
    gates: list[str],
) -> dict[str, Any]:
    family_ok = report.get("source_family") == SOURCE_FAMILY
    all_ok = report.get("all_passed") is True
    flags = {flag: report.get(flag) is True for flag in required_flags}
    passed = bool(family_ok and all_ok and all(flags.values()))
    if not passed:
        _gate(gates, f"accepted_receipt_unproved:{name}")
    return {
        "path": name,
        "sha256": RECORD_HASHES.get(name),
        "all_passed": all_ok,
        "source_family_matches": family_ok,
        "required_flags": flags,
        "passed": passed,
    }


def _provider_witness(report: Mapping[str, Any], gates: list[str]) -> dict[str, Any]:
    component = report.get("component_result") if isinstance(report.get("component_result"), Mapping) else {}
    scope = report.get("scope") if isinstance(report.get("scope"), Mapping) else {}
    flags = {
        "passed": component.get("passed") is True,
        "source_only": component.get("source_only") is True,
        "nodes_digest": component.get("nodes_digest") == EXPECTED_NODES_DIGEST,
        "node_count": component.get("node_count") == EXPECTED_NODE_COUNT,
        "ordinary_Z_orders": component.get("ordinary_Z_orders") == [0, 1, 2],
        "ordinary_Z_not_temporal_n": component.get("ordinary_Z_not_temporal_n") is True,
        "first_ten_finite_N_field_source_composition_proved": component.get(
            "first_ten_finite_N_field_source_composition_proved"
        ) is True,
        "first_ten_density_C2_NC5_identity_proved": component.get(
            "first_ten_density_C2_NC5_identity_proved"
        ) is True,
        "generic_density_identity_reused": component.get("generic_density_identity_reused") is True,
        "generic_density_identity_proved": component.get("generic_density_identity_proved") is True,
        "fixed_source_parameter_identity_proved": component.get("qualified_source_parameter_ports", {}).get("passed") is True,
        "strong_phase_equivalence_proved": component.get("strong_phase_equivalence_theorem", {}).get("passed") is True,
        "saved_exact_FZZ_formula_proved": component.get("source_recipe_theorem", {}).get("c2_density", {}).get("FZZ") is True,
        "numeric_values_computed": component.get("numeric_values_computed") is False,
        "old_operators_replayed": component.get("old_operators_replayed") is False,
    }
    scope_flags = {
        "same_density_certificate": scope.get("same_density_certificate") is False,
        "finite_N_history_identity": scope.get("finite_N_history_identity") is False,
        "physical_Q_y_Z1": scope.get("physical_Q_y_Z1") is False,
        "global_five_moment_closure": scope.get("global_five_moment_closure") is False,
        "full_goal_complete": scope.get("full_goal_complete") is False,
    }
    passed = all(flags.values()) and all(scope_flags.values())
    if not passed:
        _gate(gates, "accepted_provider_report_unproved")
    return {
        "report": PROVIDER_REPORT,
        "sha256": PROVIDER_REPORT_SHA256,
        "reused_without_rederivation": True,
        "flags": flags,
        "scope_flags": scope_flags,
        "passed": passed,
        "strong_phase_theorem_reused": component.get("strong_phase_equivalence_theorem", {}).get("passed") is True
        if isinstance(component.get("strong_phase_equivalence_theorem"), Mapping) else False,
        "strong_phase_evidence_reused": component.get("strong_phase_equivalence_theorem", {}),
        "source_equality_offsets_reused": component.get("actual_source_equality_theorem", {}).get("chart_offsets", {})
        if isinstance(component.get("actual_source_equality_theorem"), Mapping) else {},
        "loop_B_normalization_reused": component.get("strong_phase_equivalence_theorem", {}).get("loop_B_normalization")
        if isinstance(component.get("strong_phase_equivalence_theorem"), Mapping) else None,
        "phase_origin_reused": component.get("strong_phase_equivalence_theorem", {}).get("phase_origin")
        if isinstance(component.get("strong_phase_equivalence_theorem"), Mapping) else None,
        "source_composition_evidence": component.get("source_recipe_theorem", {}),
    }


def _archive_witness(
    texts: Mapping[str, str],
    gates: list[str],
    patch_gates: list[str],
) -> dict[str, Mapping[str, Any]]:
    c2 = _json(texts.get(C2_DENSITY_ARCHIVE), gates, C2_DENSITY_ARCHIVE)
    control = _json(texts.get(CONTROL_ARCHIVE), gates, CONTROL_ARCHIVE)
    c5 = _json(texts.get(C5_ARCHIVE), gates, C5_ARCHIVE)
    patch = _json(texts.get(PATCH_ARCHIVE), gates, PATCH_ARCHIVE)
    receipts = {
        C2_DENSITY_CHECK: _check_receipt(
            _json(texts.get(C2_DENSITY_CHECK), gates, C2_DENSITY_CHECK),
            (
                "current_original_17_chart_C2_signed_density_functions_installed",
                "current_original_17_chart_C2_own_rate_transport_and_five_target_functions_installed",
            ), C2_DENSITY_CHECK, gates,
        ),
        C5_CHECK: _check_receipt(
            _json(texts.get(C5_CHECK), gates, C5_CHECK),
            (
                "current_original_true_prescribed_ten_chart_signed_loop_C5_functions_installed",
                "current_original_true_prescribed_ten_chart_finite_N_density_C5_functions_installed",
                "current_original_true_prescribed_rminus_to_Rm_correction_history_C5_functions_installed",
            ), C5_CHECK, gates,
        ),
        BRIDGE_CHECK: _check_receipt(
            _json(texts.get(BRIDGE_CHECK), gates, BRIDGE_CHECK),
            (
                "current_original_true_prescribed_micro_macro_C6_functions_installed",
                "current_original_true_prescribed_bridge_Q_inertial_p2_C5_functions_installed",
            ), BRIDGE_CHECK, gates,
        ),
        R100_CHECK: _check_receipt(
            _json(texts.get(R100_CHECK), gates, R100_CHECK),
            (
                "current_original_true_prescribed_R100_R110_C6_functions_installed",
                "current_original_true_prescribed_finite_long_reference_restore_Rm_C6_functions_installed",
            ), R100_CHECK, gates,
        ),
        PATCH_CHECK: _check_receipt(
            _json(texts.get(PATCH_CHECK), gates, PATCH_CHECK),
            (
                "current_original_true_Rm_patch_finite_N_loop_density_history_C5_functions_installed",
                "current_original_true_Rm_patch_finite_N_C5_transport_installed",
            ), PATCH_CHECK, patch_gates,
        ),
    }

    expected_archive_shapes = {
        C2_DENSITY_ARCHIVE: (
            c2,
            ("exact_graph_nodes", "actual_C2_density_transport_and_target_functions"),
            12117,
        ),
        CONTROL_ARCHIVE: (
            control,
            ("exact_function_graph_nodes", "exact_current_integral_and_control_graph"),
            5607,
        ),
        C5_ARCHIVE: (
            c5,
            ("exact_graph_nodes", "prescribed_loop_C5_bindings"),
            121834,
        ),
        PATCH_ARCHIVE: (
            patch,
            ("exact_graph_nodes", "Rm_patch_finite_N_C5_bindings"),
            205147,
        ),
    }
    shape_rows: dict[str, Any] = {}
    for name, (value, keys, expected_nodes) in expected_archive_shapes.items():
        graph = value.get(keys[0])
        binding = value.get(keys[1])
        graph_ok = isinstance(graph, list) and len(graph) == expected_nodes
        binding_ok = isinstance(binding, Mapping)
        family_ok = value.get("source_family") == SOURCE_FAMILY
        passed = bool(graph_ok and binding_ok and family_ok)
        if not passed:
            _gate(patch_gates if name == PATCH_ARCHIVE else gates, f"accepted_archive_shape_unproved:{name}")
        shape_rows[name] = {
            "path": name,
            "sha256": RECORD_HASHES.get(name),
            "graph_node_count": len(graph) if isinstance(graph, list) else None,
            "expected_graph_node_count": expected_nodes,
            "binding_keys": list(binding.keys()) if isinstance(binding, Mapping) else [],
            "source_family_matches": family_ok,
            "passed": passed,
        }
    return {
        "c2": c2,
        "control": control,
        "c5": c5,
        "patch": patch,
        "receipts": receipts,
        "shape_rows": shape_rows,
    }


def _cell_token(value: Any) -> Any:
    if isinstance(value, list):
        return tuple(value)
    if isinstance(value, tuple):
        return value
    return value


def _partition_audit(
    c2_rows: Mapping[str, Mapping[str, Any]],
    control_rows: Mapping[str, Mapping[str, Any]],
    gates: list[str],
) -> dict[str, Any]:
    rows: dict[str, Any] = {}
    previous_right: Any = None
    for chart in FIRST_TEN_CHARTS:
        c2 = c2_rows.get(chart, {})
        control = control_rows.get(chart, {})
        c2_cells = c2.get("original_physical_measure_and_partitions")
        control_cells = control.get("actual_current_source_cells")
        expected = PARTITIONS.get(chart, ())
        expected_cells = [
            (expected[index], expected[index + 1])
            for index in range(max(0, len(expected) - 1))
        ]
        c2_pairs = [
            (_cell_token(item.get("exact_left")), _cell_token(item.get("exact_right")))
            for item in c2_cells
        ] if isinstance(c2_cells, list) else []
        control_pairs = [
            (_cell_token(item.get("exact_left")), _cell_token(item.get("exact_right")))
            for item in control_cells
        ] if isinstance(control_cells, list) else []
        shared_cells = c2_pairs == control_pairs == expected_cells
        adjacency = False
        if isinstance(control_cells, list) and control_cells:
            adjacency = (
                control_cells[0].get("left_offset") == control.get("left_offset")
                and control_cells[-1].get("right_offset") == control.get("right_offset")
                and all(
                    control_cells[index].get("right_offset")
                    == control_cells[index + 1].get("left_offset")
                    for index in range(len(control_cells) - 1)
                )
            )
        cross_chart_seam = previous_right is None or previous_right == control.get("left_offset")
        row_ok = bool(
            c2.get("chart") == chart
            and control.get("chart") == chart
            and shared_cells
            and adjacency
            and cross_chart_seam
            and c2.get("predecessor_ZZ_memory_retained") is True
            and c2.get("quiet_local_zero_does_not_reset_incoming") is True
            and control.get("predecessor_memory_never_reset") is True
            and control.get("quiet_local_source_zero") is False
        )
        if not row_ok:
            _gate(gates, f"first_ten_partition_coverage_unproved:{chart}")
        rows[chart] = {
            "chart": chart,
            "expected_partition": [list(pair) for pair in expected_cells],
            "C2_cells": c2_cells if isinstance(c2_cells, list) else [],
            "control_cells": control_cells if isinstance(control_cells, list) else [],
            "same_C2_and_control_cells": c2_pairs == control_pairs,
            "same_as_prescribed_native_partition": shared_cells,
            "offset_adjacency_and_telescope": adjacency,
            "previous_chart_right_is_current_left": cross_chart_seam,
            "passed": row_ok,
        }
        previous_right = control.get("right_offset")
    return rows


def _graph_children(node: Mapping[str, Any]) -> list[int]:
    children: list[int] = []
    for key, value in node.items():
        if key == "source_family":
            continue
        if isinstance(value, int) and key not in {"numerator", "denominator"}:
            # Node references use integer fields, but constants/rational
            # numerators and denominators are not graph references.
            if key in {"argument", "left", "right", "expression", "variable", "value"}:
                children.append(value)
        elif isinstance(value, list) and all(isinstance(item, int) for item in value):
            children.extend(value)
    return children


def _graph_has(node_index: Any, graph: Sequence[Any], predicate: Any) -> bool:
    if not isinstance(node_index, int) or not isinstance(graph, Sequence) or not (0 <= node_index < len(graph)):
        return False
    seen: set[int] = set()
    stack = [node_index]
    while stack:
        index = stack.pop()
        if index in seen or not (0 <= index < len(graph)):
            continue
        seen.add(index)
        value = graph[index]
        if isinstance(value, Mapping) and predicate(value):
            return True
        if isinstance(value, Mapping):
            stack.extend(_graph_children(value))
    return False


def _product_signature(graph: Sequence[Any], index: Any, *, depth: int = 0) -> tuple[Fraction, tuple[int, ...]] | None:
    """Normalize only symbolic products and exact rational coefficients."""
    if (not isinstance(index, int) or isinstance(index, bool) or not 0 <= index < len(graph)
            or depth > 128 or not isinstance(graph[index], Mapping)):
        return None
    node = graph[index]
    operation = node.get("operation")
    if operation == "exact_rational":
        try:
            return Fraction(node["numerator"], node["denominator"]), ()
        except (KeyError, TypeError, ValueError, ZeroDivisionError):
            return None
    if operation == "negative":
        value = _product_signature(graph, node.get("argument"), depth=depth+1)
        return (-value[0], value[1]) if value is not None else None
    if operation == "product":
        args = node.get("arguments")
        if not isinstance(args, list):
            return None
        coefficient, atoms = Fraction(1), []
        for arg in args:
            value = _product_signature(graph, arg, depth=depth+1)
            if value is None:
                return None
            coefficient *= value[0]
            atoms.extend(value[1])
        return coefficient, tuple(sorted(atoms))
    return Fraction(1), (index,)


def _difference_matches(graph: Sequence[Any], index: Any, left: Any, right: Any) -> bool:
    right_form = _product_signature(graph, right)
    if right_form == (Fraction(0), ()):
        return index == left
    if not isinstance(index, int) or not 0 <= index < len(graph) or not isinstance(graph[index], Mapping):
        return False
    node = graph[index]
    args = node.get("arguments")
    return bool(node.get("operation") == "sum" and isinstance(args, list) and len(args) == 2
        and args[0] == left and right_form is not None
        and _product_signature(graph, args[1]) == (-right_form[0], right_form[1]))


def _memory_rate_audit(control: Mapping[str, Any], control_rows: Mapping[str, Mapping[str, Any]],
        c5: Mapping[str, Any], gates: list[str]) -> dict[str, Any]:
    """Check actual memory exponents; no graph function is evaluated."""
    c2_graph = control.get("exact_function_graph_nodes", [])
    c5_graph = c5.get("exact_graph_nodes", [])
    charts = c5.get("prescribed_loop_C5_bindings", {}).get("charts", {})
    out: dict[str, Any] = {}
    for chart in FIRST_TEN_CHARTS:
        old = control_rows.get(chart, {})
        current = charts.get(chart, {})
        width_form = _product_signature(c2_graph, old.get("width"))
        checks: dict[str, Any] = {}
        for field, rate in RATE_VALUES.items():
            index = old.get("memory", {}).get(field)
            c2_ok = _product_signature(c2_graph, index) == (Fraction(1), ()) if rate == 0 else False
            if rate and isinstance(index, int) and 0 <= index < len(c2_graph):
                node = c2_graph[index]
                c2_ok = (isinstance(node, Mapping) and node.get("operation") == "exp"
                    and width_form is not None and _product_signature(c2_graph, node.get("argument"))
                    == (-rate*width_form[0], width_form[1]))
            index5 = current.get("memory", {}).get(field)
            c5_ok = _product_signature(c5_graph, index5) == (Fraction(1), ()) if rate == 0 else False
            if rate and isinstance(index5, int) and 0 <= index5 < len(c5_graph):
                node5 = c5_graph[index5]
                exponent = _product_signature(c5_graph, node5.get("argument")) if isinstance(node5, Mapping) else None
                measure = _product_signature(c5_graph, current.get("physical_measure"))
                domain = current.get("domain")
                candidates = list(exponent[1]) if exponent is not None else []
                if isinstance(domain, list) and len(domain) == 2:
                    candidates.append(current.get("coordinate"))
                    for delta in candidates:
                        if not _difference_matches(c5_graph, delta, current.get("coordinate"), domain[0]):
                            continue
                        delta_form = _product_signature(c5_graph, delta)
                        if measure is not None and delta_form is not None:
                            expected = (-rate*measure[0]*delta_form[0], tuple(sorted(measure[1]+delta_form[1])))
                            c5_ok = bool(node5.get("operation") == "exp" and exponent == expected)
                            if c5_ok:
                                break
            checks[field] = dict(rate=str(rate), C2_actual_exponent_matches=c2_ok,
                C5_actual_exponent_matches=c5_ok, passed=bool(c2_ok and c5_ok))
        passed = all(row["passed"] for row in checks.values())
        if not passed:
            _gate(gates, "actual_own_rate_memory_unproved:"+chart)
        out[chart] = dict(checks=checks, passed=passed,
            inspection="Symbolic exact rational rate coefficients and opaque function-port products only; no field values computed")
    return out


def _jacobian_audit(
    control: Mapping[str, Any],
    control_rows: Mapping[str, Mapping[str, Any]],
    gates: list[str],
    *,
    source_geometry_proved: bool,
) -> dict[str, Any]:
    graph = control.get("exact_function_graph_nodes")
    expected = {
        "first_micro": {"parameter": "hbB"},
        "second_micro": {"parameter": "hbB"},
        "frozen_macro": {"parameters": ("hbB",), "contains_constant": -2},
        "first_switch": {"parameter": "hbS"},
        "second_switch": {"parameter": "hbS"},
        "post_power": {"parameters": ("hbS",), "contains_constant": -2},
        "long_reshape": {"parameter": "T"},
        "reference": {"parameters": ("T", "logC"), "contains_constant": -8},
        "restoration": {"rational": (1, 1)},
        "postrestore": {"rational": (1, 1)},
    }
    rows: dict[str, Any] = {}
    for chart in FIRST_TEN_CHARTS:
        row = control_rows.get(chart, {})
        index = row.get("Jacobian")
        node = graph[index] if isinstance(graph, Sequence) and isinstance(index, int) and 0 <= index < len(graph) else None
        spec = expected[chart]
        parameter = spec.get("parameter")
        parameter_ok = bool(
            parameter is None
            or isinstance(node, Mapping)
            and node.get("operation") == "current_original_source_parameter"
            and node.get("name") == parameter
        )
        if "parameters" in spec:
            parameter_ok = all(
                _graph_has(index, graph, lambda value, name=name: value.get("name") == name)
                for name in spec["parameters"]
            ) if isinstance(graph, Sequence) else False
        constant_ok = True
        if "contains_constant" in spec:
            target = spec["contains_constant"]
            constant_ok = _graph_has(
                index, graph,
                lambda value: value.get("numerator") == target and value.get("denominator") == 1,
            ) if isinstance(graph, Sequence) else False
        rational_ok = True
        if "rational" in spec:
            rational_ok = isinstance(node, Mapping) and node.get("operation") == "exact_rational" and (
                node.get("numerator"), node.get("denominator")
            ) == spec["rational"]
        row_ok = bool(
            isinstance(index, int)
            and isinstance(node, Mapping)
            and source_geometry_proved
            and parameter_ok
            and constant_ok
            and rational_ok
        )
        if not row_ok:
            _gate(gates, f"physical_jacobian_source_unproved:{chart}")
        rows[chart] = {
            "Jacobian_node": index,
            "node_operation": node.get("operation") if isinstance(node, Mapping) else None,
            "source_parameter_or_rational": (
                node.get("name") if isinstance(node, Mapping) else None
            ),
            "expected_formula": GEOMETRY[chart]["Jacobian"],
            "source_exact_geometry_and_sy_diff_proved": source_geometry_proved,
            "passed": row_ok,
        }
    return rows


def _history_audit(
    c2_rows: Mapping[str, Mapping[str, Any]],
    control_rows: Mapping[str, Mapping[str, Any]],
    c5_bindings: Mapping[str, Any],
    gates: list[str],
) -> tuple[dict[str, Any], bool]:
    rows: dict[str, Any] = {}
    previous_c2_outgoing: Mapping[str, Any] | None = None
    previous_control_outgoing: Mapping[str, Any] | None = None
    for chart in FIRST_TEN_CHARTS:
        c2 = c2_rows.get(chart, {})
        control = control_rows.get(chart, {})
        checks: dict[str, bool] = {}
        checks["all_rate_ports_present"] = all(
            isinstance(c2.get(key), Mapping) and set(c2[key]) == set(RATES)
            for key in ("density", "incoming", "contributions", "outgoing", "original_own_rate_memory")
        ) and all(
            isinstance(control.get(key), Mapping) and set(control[key]) == set(RATES)
            for key in ("density", "incoming", "contributions", "memory", "outgoing")
        )
        checks["C2_old_prefix_binds_control_C1_rows"] = checks["all_rate_ports_present"] and all(
            isinstance(c2["density"].get(rate), list)
            and isinstance(control["density"].get(rate), list)
            and c2["density"][rate][:2] == control["density"][rate]
            and isinstance(c2["incoming"].get(rate), list)
            and c2["incoming"][rate][:2] == control["incoming"][rate]
            and isinstance(c2["contributions"].get(rate), list)
            and c2["contributions"][rate][:2] == control["contributions"][rate]
            and isinstance(c2["outgoing"].get(rate), list)
            and c2["outgoing"][rate][:2] == control["outgoing"][rate]
            for rate in RATES
        ) if checks["all_rate_ports_present"] else False
        checks["C2_has_exact_Z012_rows"] = checks["all_rate_ports_present"] and all(
            len(c2["density"][rate]) == 3
            and len(c2["incoming"][rate]) == 3
            and len(c2["contributions"][rate]) == 3
            and len(c2["outgoing"][rate]) == 3
            for rate in RATES
        ) if checks["all_rate_ports_present"] else False
        checks["C2_memory_equals_control_memory"] = checks["all_rate_ports_present"] and all(
            c2["original_own_rate_memory"][rate] == control["memory"][rate]
            for rate in RATES
        ) if checks["all_rate_ports_present"] else False
        checks["rate_zero_pressure_memory_is_one"] = (
            c2.get("original_own_rate_memory", {}).get("p") == 1
            and control.get("memory", {}).get("p") == 1
        )
        if previous_c2_outgoing is None:
            checks["first_incoming_is_zero_source_port"] = all(
                isinstance(c2.get("incoming", {}).get(rate), list)
                and c2["incoming"][rate][:2] == [0, 0]
                and c2["incoming"][rate][2] == 0
                for rate in RATES
            )
        else:
            checks["incoming_is_previous_C2_outgoing"] = all(
                c2.get("incoming", {}).get(rate) == previous_c2_outgoing.get(rate)
                for rate in RATES
            )
        if previous_control_outgoing is None:
            checks["first_control_incoming_is_zero_source_port"] = all(
                control.get("incoming", {}).get(rate) == [0, 0] for rate in RATES
            )
        else:
            checks["control_incoming_is_previous_outgoing"] = all(
                control.get("incoming", {}).get(rate) == previous_control_outgoing.get(rate)
                for rate in RATES
            )
        checks["control_phase_is_Z_independent"] = control.get("loop", {}).get("phase_Z_exact_zero") is True
        checks["C2_quiet_branch_does_not_reset"] = (
            c2.get("quiet_local_zero_does_not_reset_incoming") is True
            and control.get("predecessor_memory_never_reset") is True
        )
        previous_c2_outgoing = c2.get("outgoing") if isinstance(c2.get("outgoing"), Mapping) else None
        previous_control_outgoing = control.get("outgoing") if isinstance(control.get("outgoing"), Mapping) else None

        c5 = c5_bindings.get(chart) if isinstance(c5_bindings, Mapping) else None
        endpoint = c5.get("history") if isinstance(c5, Mapping) else None
        c5_endpoint = c5_bindings.get("__endpoints__", {}).get(chart) if isinstance(c5_bindings.get("__endpoints__"), Mapping) else None
        checks["C5_history_port_present"] = isinstance(endpoint, Mapping) and set(endpoint) == set(RATES)
        checks["C5_endpoint_substitution_port_present"] = isinstance(c5_endpoint, Mapping) and set(c5_endpoint) == set(RATES)
        checks["chart_passed"] = all(checks.values())
        if not checks["chart_passed"]:
            _gate(gates, f"first_ten_nonreset_history_unproved:{chart}")
        rows[chart] = {
            "chart": chart,
            "checks": checks,
            "C2_density_Z012": c2.get("density", {}),
            "C2_incoming_Z012": c2.get("incoming", {}),
            "C2_contributions_Z012": c2.get("contributions", {}),
            "C2_outgoing_Z012": c2.get("outgoing", {}),
            "C2_original_own_rate_memory": c2.get("original_own_rate_memory", {}),
            "control_old_C1_density": control.get("density", {}),
            "control_old_C1_incoming": control.get("incoming", {}),
            "control_old_C1_contributions": control.get("contributions", {}),
            "control_old_C1_memory": control.get("memory", {}),
            "control_old_C1_outgoing": control.get("outgoing", {}),
            "C5_history_source_ports": endpoint,
            "C5_endpoint_ports": c5_endpoint,
        }
    return rows, bool(rows) and all(row["checks"]["chart_passed"] for row in rows.values())


def _interior_specs(
    history_rows: Mapping[str, Any],
    partition_rows: Mapping[str, Any],
    control_rows: Mapping[str, Mapping[str, Any]],
    c2_rows: Mapping[str, Mapping[str, Any]],
    gates: list[str],
) -> dict[str, Any]:
    """Expose a source-defined interior Duhamel specification per chart."""

    result: dict[str, Any] = {}
    for chart in FIRST_TEN_CHARTS:
        c2 = c2_rows.get(chart, {})
        control = control_rows.get(chart, {})
        parts = partition_rows.get(chart, {})
        row = history_rows.get(chart, {})
        cell_rows = parts.get("control_cells", [])
        density = c2.get("density", {})
        source_cells = []
        for cell in cell_rows:
            source_cells.append({
                "exact_left": cell.get("exact_left"),
                "exact_right": cell.get("exact_right"),
                "lower": cell.get("lower"),
                "upper": cell.get("upper"),
                "left_offset": cell.get("left_offset"),
                "right_offset": cell.get("right_offset"),
                "physical_measure": "control.Jacobian * dt",
                "Jacobian_applied_once": True,
            })
        passed = bool(
            row.get("checks", {}).get("chart_passed") is True
            and parts.get("passed") is True
            and isinstance(density, Mapping) and set(density) == set(RATES)
            and isinstance(control.get("Jacobian"), int)
            and isinstance(control.get("offset"), int)
            and isinstance(control.get("right_offset"), int)
        )
        if not passed:
            _gate(gates, f"first_ten_interior_source_unproved:{chart}")
        result[chart] = {
            "chart": chart,
            "source_function": "New continuous C2 Duhamel definition from the actual density Z012 and incoming ports; saved C2 history ports are endpoints",
            "density_Z012_ports": density,
            "incoming_Z012_ports": c2.get("incoming", {}),
            "contribution_Z012_ports": c2.get("contributions", {}),
            "outgoing_Z012_ports": c2.get("outgoing", {}),
            "own_rate_memory": c2.get("original_own_rate_memory", {}),
            "transport_formula": (
                "H_out_Zq=exp(-rate*(right_offset-left_offset))*H_in_Zq + "
                "sum_native_partition definite_integral(" 
                "exp(-rate*(right_offset-offset(t)))*density_Zq(t)*Jacobian dt), q=0,1,2"
            ),
            "continuous_interior_definition": {
                "operation": "source_function_definition",
                "new_definition_not_an_archived_interior_H_port": True,
                "formula": "H_Zq(s)=exp(-rate*(offset(s)-left_offset))*H_in_Zq + integral_lower^s Jacobian*exp(-rate*(offset(s)-offset(t)))*density_Zq(t) dt",
                "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
                "density_archive": C2_DENSITY_ARCHIVE,
                "geometry_archive": CONTROL_ARCHIVE,
                "density_Z012": density,
                "incoming_Z012": c2.get("incoming", {}),
                "offset_node": control.get("offset"),
                "left_offset_node": control.get("left_offset"),
                "Jacobian_node": control.get("Jacobian"),
                "domain_nodes": control.get("domain"),
                "native_variable": "native_" + chart,
                "endpoint_identification_basis": "Finite integral additivity over the complete native partition, with the same incoming memory and right offset",
            },
            "rate_zero_memory": "memory=1 and incoming is retained",
            "source_cells": source_cells,
            "chart_domain": GEOMETRY[chart]["domain"],
            "affine_offset": GEOMETRY[chart]["offset"],
            "physical_Jacobian": GEOMETRY[chart]["Jacobian"],
            "phase_is_Z_independent": control.get("loop", {}).get("phase_Z_exact_zero") is True,
            "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
            "endpoint_derivative_source_spec": {
                "operator": "D_y=Jacobian^{-1}*partial_native_coordinate",
                "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
                "formula": "D_y H_Zq(right)=density_Zq(right)-rate*H_out_Zq(right)",
                "right_endpoint_density_source": "function_substitution(density_Zq, coordinate, upper)",
                "uses_right_endpoint_outgoing_port": True,
                "incoming_used_only_at_left_boundary": True,
                "patch_inlet_derivative_not_claimed": True,
            },
            "passed": passed,
        }
    return result


def _endpoint_audit(
    c5_binding: Mapping[str, Any],
    c5_graph: Sequence[Any],
    c2_rows: Mapping[str, Mapping[str, Any]],
    gates: list[str],
) -> dict[str, Any]:
    endpoints = c5_binding.get("endpoints") if isinstance(c5_binding, Mapping) else {}
    charts = c5_binding.get("charts") if isinstance(c5_binding, Mapping) else {}
    rows: dict[str, Any] = {}
    previous_outgoing: Mapping[str, Any] | None = None
    previous_c5_endpoint: Mapping[str, Any] | None = None
    for chart in FIRST_TEN_CHARTS:
        chart_row = charts.get(chart) if isinstance(charts, Mapping) else None
        endpoint = endpoints.get(chart) if isinstance(endpoints, Mapping) else None
        outgoing = c2_rows.get(chart, {}).get("outgoing")
        current_ok = isinstance(chart_row, Mapping) and isinstance(endpoint, Mapping) and set(endpoint) == set(RATES)
        chain_ok = previous_outgoing is None or all(
            c2_rows.get(chart, {}).get("incoming", {}).get(rate) == previous_outgoing.get(rate)
            for rate in RATES
        )
        incoming = chart_row.get("incoming", {}) if isinstance(chart_row, Mapping) else {}
        history = chart_row.get("history", {}) if isinstance(chart_row, Mapping) else {}
        coordinate = chart_row.get("coordinate") if isinstance(chart_row, Mapping) else None
        domain = chart_row.get("domain") if isinstance(chart_row, Mapping) else None
        if previous_c5_endpoint is None:
            c5_chain_ok = (isinstance(incoming, Mapping) and set(incoming) == set(RATES)
                and all(isinstance(incoming.get(rate), list) and incoming[rate] == [0]*6 for rate in RATES)
                and isinstance(c5_graph, Sequence) and len(c5_graph) > 0
                and isinstance(c5_graph[0], Mapping)
                and c5_graph[0].get("operation") == "exact_rational"
                and c5_graph[0].get("numerator") == 0)
        else:
            c5_chain_ok = incoming == previous_c5_endpoint
        substitutions: dict[str, bool] = {}
        for rate in RATES:
            qports = endpoint.get(rate) if isinstance(endpoint, Mapping) else None
            hports = history.get(rate) if isinstance(history, Mapping) else None
            shape_ok = (isinstance(qports, list) and len(qports) == 6
                and isinstance(hports, list) and len(hports) == 6
                and isinstance(domain, list) and len(domain) == 2)
            substitutions[rate] = bool(shape_ok and all(
                isinstance(index, int) and not isinstance(index, bool) and 0 <= index < len(c5_graph)
                and isinstance(c5_graph[index], Mapping)
                and c5_graph[index].get("operation") == "function_substitution"
                and c5_graph[index].get("expression") == hports[q]
                and c5_graph[index].get("variable") == coordinate
                and c5_graph[index].get("value") == domain[1]
                and c5_graph[index].get("Z_independent_substitution") is True
                for q, index in enumerate(qports)
            ))
        passed = bool(current_ok and chain_ok and isinstance(outgoing, Mapping)
            and c5_chain_ok and all(substitutions.values()))
        if not passed:
            _gate(gates, f"endpoint_function_port_chain_unproved:{chart}")
        rows[chart] = {
            "chart": chart,
            "endpoint_ports": endpoint,
            "endpoint_rows_are_values": False,
            "endpoint_operation": "a.at(history,coordinate,upper)",
            "C2_incoming_is_previous_outgoing": chain_ok,
            "C5_incoming_is_previous_endpoint_or_original_zero": c5_chain_ok,
            "C5_endpoint_substitution_nodes_match_all_six_orders": substitutions,
            "becomes_next_incoming": chain_ok and c5_chain_ok,
            "C2_outgoing_port": outgoing,
            "passed": passed,
        }
        previous_outgoing = outgoing if isinstance(outgoing, Mapping) else None
        previous_c5_endpoint = endpoint if isinstance(endpoint, Mapping) else None
    return rows


def _patch_witness(
    patch: Mapping[str, Any],
    c5_binding: Mapping[str, Any],
    c2_rows: Mapping[str, Mapping[str, Any]],
    gates: list[str],
) -> dict[str, Any]:
    binding = patch.get("Rm_patch_finite_N_C5_bindings") if isinstance(patch, Mapping) else {}
    chart = binding.get("chart") if isinstance(binding, Mapping) else {}
    definition = binding.get("definition") if isinstance(binding, Mapping) else {}
    incoming = binding.get("incoming") if isinstance(binding, Mapping) else {}
    postrestore = c2_rows.get("postrestore", {}).get("outgoing")
    checks = {
        "patch_binding_present": isinstance(binding, Mapping),
        "patch_chart_present": isinstance(chart, Mapping),
        "patch_nonzero_incoming_port_present": isinstance(incoming, Mapping) and set(incoming) == set(RATES),
        "patch_definition_retains_nonzero_Rm_inlet": definition.get(
            "nonzero_true_Rm_correction_and_signed_pressure_inlets_retained"
        ) is True,
        "patch_definition_retains_log_transport": definition.get(
            "exact_log_coordinate_transport_and_rate_zero_pressure_memory"
        ) is True,
        "postrestore_C2_outgoing_retained": isinstance(postrestore, Mapping) and set(postrestore) == set(RATES),
        "incoming_is_saved_C5_postrestore_endpoint": incoming == c5_binding.get("Rm_correction")
            == c5_binding.get("endpoints", {}).get("postrestore"),
    }
    if not all(checks.values()):
        _gate(gates, "patch_transport_port_unproved")
    # The finite-N patch source is admitted as a retained incoming-port recipe;
    # its independent patch-source identity theorem remains intentionally open.
    return {
        "checks": checks,
        "passed": all(checks.values()),
        "patch_source_theorem_proved": False,
        "patch_identity_stays_false": True,
        "patch_incoming_source": "C5 prescribed_loop_C5_functions.Rm_correction",
        "postrestore_outgoing_source": "C2_density_transport.windows[postrestore].outgoing",
        "postrestore_outgoing_Z012": postrestore,
        "patch_receipt_unresolved_gates_are_separate": True,
    }


def _module_hash() -> str | None:
    # The runner supplies the implementation hash separately.  Avoid reading
    # the module from inside derive, which keeps this consumer source-only.
    return None


def derive(
    nodes: Sequence[Any],
    source_records: Mapping[str, Any],
    *,
    expected_nodes_digest: str,
) -> dict[str, Any]:
    """Inspect authenticated source records and bind nonreset first-ten history."""

    gates: list[str] = []
    patch_gates: list[str] = []
    actual_digest = ""
    try:
        actual_digest = canonical_digest(nodes)
    except (TypeError, ValueError, OverflowError):
        _gate(gates, "nodes_digest_error")
    valid_nodes = isinstance(nodes, Sequence) and not isinstance(nodes, (str, bytes, bytearray))
    if expected_nodes_digest != EXPECTED_NODES_DIGEST:
        _gate(gates, "expected_nodes_digest_mismatch")
    if actual_digest != EXPECTED_NODES_DIGEST:
        _gate(gates, "nodes_digest_mismatch")
    if not valid_nodes or len(nodes) != EXPECTED_NODE_COUNT:
        _gate(gates, "node_count_mismatch")
    if not isinstance(source_records, Mapping):
        _gate(gates, "source_records_not_mapping")
        records: Mapping[str, Any] = {}
    else:
        records = source_records

    texts: dict[str, str] = {}
    for name in RECORD_HASHES:
        if name.endswith(".gz"):
            text = _raw_archive(records, name, gates)
        else:
            text = _record(records, name, gates)
        if text is not None:
            texts[name] = text

    provider_report = _json(texts.get(PROVIDER_REPORT), gates, PROVIDER_REPORT)
    provider = _provider_witness(provider_report, gates)

    # Every source file is parsed as plain text only.  No imported module or
    # helper is loaded, and no function body is invoked.
    source_names = (
        C2_DENSITY, CONTROL_CONNECTION, C5_SOURCE, CORE_INLET,
        BRIDGE_SOURCE, R100_SOURCE, PATCH_SOURCE,
    )
    trees = {name: _tree(texts.get(name), name, gates) for name in source_names}
    source_ast = _source_ast_witness(trees, gates, patch_gates)

    archives = _archive_witness(texts, gates, patch_gates)
    c2_archive = archives["c2"]
    control_archive = archives["control"]
    c5_archive = archives["c5"]
    patch_archive = archives["patch"]
    c2_binding = c2_archive.get("actual_C2_density_transport_and_target_functions", {})
    control_binding = control_archive.get("exact_current_integral_and_control_graph", {})
    c2_list = c2_binding.get("windows") if isinstance(c2_binding, Mapping) else None
    control_list = control_binding.get("windows") if isinstance(control_binding, Mapping) else None
    c5_binding = c5_archive.get("prescribed_loop_C5_bindings", {})
    c5_charts = c5_binding.get("charts") if isinstance(c5_binding, Mapping) else {}

    c2_rows = {
        row.get("chart"): row for row in c2_list
        if isinstance(row, Mapping) and isinstance(row.get("chart"), str)
    } if isinstance(c2_list, list) else {}
    control_rows = {
        row.get("chart"): row for row in control_list
        if isinstance(row, Mapping) and isinstance(row.get("chart"), str)
    } if isinstance(control_list, list) else {}

    chart_order_ok = (
        isinstance(c2_list, list)
        and [row.get("chart") for row in c2_list[:10] if isinstance(row, Mapping)] == list(FIRST_TEN_CHARTS)
        and isinstance(control_list, list)
        and [row.get("chart") for row in control_list[:10] if isinstance(row, Mapping)] == list(FIRST_TEN_CHARTS)
        and isinstance(c5_charts, Mapping)
        and tuple(name for name in c5_charts if name in FIRST_TEN_CHARTS) == FIRST_TEN_CHARTS
    )
    if not chart_order_ok:
        _gate(gates, "first_ten_chart_order_unproved")

    partitions = _partition_audit(c2_rows, control_rows, gates)
    geometry_ast = source_ast.get("control_geometry_and_graph", {})
    source_geometry_proved = bool(
        geometry_ast.get("exact_geometry_function_present") is True
        and geometry_ast.get("offsets_all_match") is True
        and geometry_ast.get("control_graph_function_present") is True
        and geometry_ast.get("Jacobian_assignment_is_symbolic_offset_derivative") is True
    )
    jacobians = _jacobian_audit(
        control_archive,
        control_rows,
        gates,
        source_geometry_proved=source_geometry_proved,
    )
    # Make endpoint references available beside the C5 chart rows without
    # mutating source input data.
    c5_for_history = dict(c5_charts) if isinstance(c5_charts, Mapping) else {}
    c5_for_history["__endpoints__"] = c5_binding.get("endpoints", {}) if isinstance(c5_binding, Mapping) else {}
    histories, histories_passed = _history_audit(c2_rows, control_rows, c5_for_history, gates)
    memories = _memory_rate_audit(control_archive, control_rows, c5_archive, gates)
    interiors = _interior_specs(histories, partitions, control_rows, c2_rows, gates)
    endpoints = _endpoint_audit(c5_binding, c5_archive.get("exact_graph_nodes", []), c2_rows, gates)
    patch = _patch_witness(patch_archive, c5_binding, c2_rows, patch_gates)

    all_first_ten_gates_clear = (
        not gates
        and provider.get("passed") is True
        and chart_order_ok
        and histories_passed
        and all(row.get("passed") is True for row in partitions.values())
        and all(row.get("passed") is True for row in jacobians.values())
        and all(row.get("passed") is True for row in interiors.values())
        and all(row.get("passed") is True for row in endpoints.values())
        and all(row.get("passed") is True for row in memories.values())
    )
    if not all_first_ten_gates_clear:
        _gate(gates, "first_ten_nonreset_history_identity_unproved")
    endpoint_identity = bool(all_first_ten_gates_clear)

    # The source-defined physical geometry is retained as conditional source
    # reasoning.  No numerical positivity or field value is admitted here.
    geometry = {
        "charts": {
            chart: {
                **GEOMETRY[chart],
                "domain": list(GEOMETRY[chart]["domain"]),
                "native_partition": [list(pair) for pair in zip(PARTITIONS[chart], PARTITIONS[chart][1:])],
                "Jacobian_is_source_derivative": True,
                "rate_kernel_uses_same_offset": True,
            }
            for chart in FIRST_TEN_CHARTS
        },
        "positive_width_conditions": dict(POSITIVE_WIDTH_CONDITIONS),
        "positive_widths_numerically_admitted": False,
        "frozen_macro_derivative_corrected": "Y-2*hbB",
        "seams": {
            "cell_offsets_are_adjacent": all(row.get("offset_adjacency_and_telescope") for row in partitions.values()),
            "endpoint_right_offset_is_next_left_offset": all(row.get("previous_chart_right_is_current_left") for row in partitions.values()),
            "first_ten_chart_order_is_source_order": chart_order_ok,
        },
    }

    result = {
        "schema": SCHEMA,
        "source_only": True,
        "passed": bool(all_first_ten_gates_clear and not gates),
        "nodes_digest": actual_digest,
        "expected_nodes_digest": EXPECTED_NODES_DIGEST,
        "node_count": len(nodes) if valid_nodes else None,
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "ordinary_Z_not_temporal_n": True,
        "time_order_recursion_claimed": False,
        "timeorder_recursion_claimed": False,
        "coveredcharts": list(FIRST_TEN_CHARTS) if all_first_ten_gates_clear else [],
        "covered_charts": list(FIRST_TEN_CHARTS) if all_first_ten_gates_clear else [],
        "first_ten_nonreset_history_identity_proved": bool(all_first_ten_gates_clear),
        "first_ten_history_endpoint_identity_proved": endpoint_identity,
        "first_ten_density_C2_NC5_identity_proved": provider["flags"].get(
            "first_ten_density_C2_NC5_identity_proved", False,
        ),
        "first_ten_finite_N_field_source_composition_proved": provider["flags"].get(
            "first_ten_finite_N_field_source_composition_proved", False,
        ),
        "generic_density_identity_reused": provider["flags"].get("generic_density_identity_reused", False),
        "generic_density_identity_proved": provider["flags"].get("generic_density_identity_proved", False),
        "accepted_provider_report": provider,
        "source_ast_witness": source_ast,
        "accepted_source_receipts": archives["receipts"],
        "archive_shape_witness": archives["shape_rows"],
        "source_geometry_and_jacobians": geometry,
        "physical_jacobian_audit": jacobians,
        "actual_own_rate_memory_audit": memories,
        "partition_coverage": partitions,
        "first_ten_history_ports": histories,
        "first_ten_interior_source_specs": interiors,
        "perchart_history_bindings": interiors,
        "history_identity_proof": {
            "premise": "The saved actual first-ten density functions satisfy rho_C2=N*rho_C5 for ordinary Z0/1/2, with N independent of Z.",
            "base": "Both correction inlets vanish only at original Rminus, as checked by the archived zero function ports.",
            "step": "On each chart the same exact rate, affine physical logarithmic offset, and Jacobian give the same scalar Duhamel operator. Its linearity and additivity over the complete native partition carry H_C2=N*H_C5 from the incoming endpoint to the whole chart and its outgoing endpoint.",
            "induction": "Every next incoming is the previous outgoing in each archive; there is no interior reset.",
            "derivatives": "The fixed endpoints, offset, rate, N and Jacobian are independent of Z. Genuine C2 source functions on the conditional compact chart domain permit ordinary Z derivatives through the finite integrals twice.",
            "source_function_identity_only": True,
            "numerical_values_or_error_reductions_computed": False,
            "proved": bool(all_first_ten_gates_clear),
        },
        "domain_evidence": {"source_only": True, "conditional_original_common_admissible_domain": True,
            "positive_width_conditions": dict(POSITIVE_WIDTH_CONDITIONS), "ordinary_Z_orders": [0,1,2],
            "numerical_domain_admission_proved": False, "endpoint_rows_are_values": False},
        "endpoint_ports": {
            "endpoint_rows_are_values": False,
            "source_function_ports_retained": True,
            "chart_ports": endpoints,
            "a_at_endpoint_becomes_next_incoming": endpoint_identity,
        },
        "patch_provider": patch,
        "patch_unresolved_gates": sorted(set(patch_gates)
            | set(provider_report.get("component_result", {}).get("patch_unresolved_gates", []))
            | set(provider_report.get("component_result", {}).get("patch_provider", {}).get("unresolved_gates", []))
            | {"patch_source_function_identity_and_history_unproved"}),
        "patch_source_composition_proved": False,
        "patch_identity_stays_false": True,
        "numeric_values_computed": False,
        "old_operators_replayed": False,
        "mathematical_operator_replayed": False,
        "finite_N_history_identity": False,
        "same_density_certificate": False,
        "physical_Q_y_Z1": False,
        "global_five_moment_closure": False,
        "full_goal_complete": False,
        "source_inspection_audit": {
            "source_records_only": True,
            "static_source_AST_graph_and_saved_JSON_data_inspection_only": True,
            "nodes_hashed_only": True,
            "new_source_consumer_calls": 1,
            "accepted_ancestor_operator_calls": 0,
            "accepted_mathematical_operator_calls": 0,
            "old_mathematical_operator_calls": 0,
            "history_operator_calls": 0,
            "numeric_field_evaluated": False,
            "quadrature_or_norm_calls": 0,
            "recurrence_calls": 0,
            "target_calls": 0,
        },
        "accepted_source_records": {
            name: {"sha256": digest, "path": "experiments/root_st073/" + name}
            for name, digest in RECORD_HASHES.items()
        },
        "unresolved_gates": sorted(set(gates)),
        "scope_limits": {
            "history_identity_is_first_ten_only": True,
            "all_first_ten_histories_use_nonreset_incoming": histories_passed,
            "C2_ZZ_rows_are_source_ports": True,
            "ordinary_Z_orders_are_not_temporal_orders": True,
            "endpoint_rows_are_values": False,
            "physical_width_conditions_are_conditional": True,
            "patch_source_theorem_proved": False,
            "global_five_moment_closure": False,
            "full_goal_complete": False,
        },
        "record_hashes": dict(RECORD_HASHES),
        "module_sha256": _module_hash(),
        "runtime_model_metadata": None,
    }
    return result


__all__ = [
    "EXPECTED_NODE_COUNT",
    "EXPECTED_NODES_DIGEST",
    "RECORD_HASHES",
    "SCHEMA",
    "canonical_digest",
    "derive",
]
