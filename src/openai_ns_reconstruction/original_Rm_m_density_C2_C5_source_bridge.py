"""Bounded source certificate for the C2/C5 density ``m`` bridge.

This module does not import or execute a producer, a checker, a graph
builder, an integral, a recurrence, or a norm.  It authenticates the supplied
205300-node graph by its canonical digest and inspects only pinned source text
and accepted text receipts.  The admitted statement is the exact linear
source identity

``C2_m = B = N * (B_over_Pstar / N) = N * C5_m``.

The equality is admitted only when the C2 and C5 source definitions and the
same live source-provider chain are both present.  A matching field name or a
generic family name is never treated as a provider binding.  Post-restore and
patch rows are reported independently; the patch row remains unresolved when
its optional direct C5 patch source is absent.
"""

from __future__ import annotations

import ast
from collections.abc import Mapping, Sequence
import hashlib
import json
from pathlib import Path
from typing import Any


SCHEMA = "original_Rm_m_density_C2_C5_source_bridge_v1"
EXPECTED_NODES_DIGEST = "a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3"
EXPECTED_NODE_COUNT = 205300

# This is the single accepted source family.  The provider records below are
# allowed to omit this object in their compact receipts; when it is present,
# every hash must still agree exactly.
SOURCE_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}

PREFIX = "lei_ren_part1_paper_compliant_current_original_"

# Only plain UTF-8 source and receipt records are accepted.  In particular,
# no gzip archive is decoded here.  The optional patch records are pinned so a
# caller can admit the direct patch row, while their absence leaves that row
# explicitly unresolved instead of weakening the first-nine admission.
RECORD_HASHES = {
    PREFIX + "C2_phase_primitives.py": "84f5da3eb64096f82df88dd828b22b09bce18192a03925d91cfee0b193902d2c",
    PREFIX + "C2_phase_primitives_check.json": "76ed1af50ec799c49081ab67881b189b42e0d3137bcca770a67766d4111dff8b",
    PREFIX + "C2_density_transport.py": "8bfa32c87c73c5bedf3056c7171707ae934dafd950bd0df874dbd3b3450379d1",
    PREFIX + "C2_density_transport_check.json": "c4d2c62dd67e9dd7e171fabbdb9ba18308e420c5e0c98d01147c708840b8f583",
    PREFIX + "prescribed_loop_C5_functions.py": "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d",
    PREFIX + "prescribed_loop_C5_functions_check.json": "2c77e80df04fae615f8ea909fc7b606bfccfa27659cb67744ae8d84020b1a9e7",
    PREFIX + "whole_Z_bridge_source.py": "df8b204d023c1d8762da16d07190601ebd1e11447077b01a734bba8baa595a6a",
    PREFIX + "whole_Z_bridge_source_check.json": "91d9ba8119d9dba85322ea33791ac3126c22ad85dfd94bf9a961c49fc18268a3",
    PREFIX + "whole_Z_bridge_correlated_functions.py": "bb29100e671b6a1927f9897385056645a17d859108e537eec760a3904d66f832",
    PREFIX + "whole_Z_bridge_correlated_functions_check.json": "bbf33386a8a591eef946ca5f3f262dd92ac6dd1490bc8bfc7421ee1c6eca49da",
    PREFIX + "whole_Z_all_N_source_functions.py": "f0f6058bf77e52d9f5f04db3259ca9f737dc0fdb9281e2650b7a6ae5ca3f25fe",
    PREFIX + "whole_Z_all_N_source_functions_check.json": "b141175a0b476d2c7b325028ef6d45696ca5b6e52cdafcac775d2d98e2dd8ccb",
    PREFIX + "whole_Z_all_N_switch_functions.py": "8a3e9581a393fcbc59fdad7f24055b0176eb7edca8a76676a35c1b23b114922c",
    PREFIX + "whole_Z_all_N_switch_functions_check.json": "3b69567b7e51aef13be5f35492ebf8f4013b43b8d15c2d9f7d175bb4292c0497",
    PREFIX + "whole_Z_all_N_long_patch_functions.py": "847120150d0c8a83c2861f440c8f69514358b55633cfe8964f6de25bd5aa3ddc",
    PREFIX + "whole_Z_all_N_long_patch_functions_check.json": "2162a7b7fbf322a53065e9ba118e9466a7b17a3baff86d569ed35935b18e8a02",
    PREFIX + "prescribed_bridge_C6_functions.py": "8ab2c26092db06453ed5a6047d8905e37fdad17eaf9c4947048498551425c605",
    PREFIX + "prescribed_bridge_C6_functions_check.json": "9a0969a8cd4a5989dbd6272d6dd223533400113824f9352b4e57a83196bff61e",
    PREFIX + "prescribed_R100_Rm_C6_functions.py": "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32",
    PREFIX + "prescribed_R100_Rm_C6_functions_check.json": "494caf482bffe5475c011a35a1b7d95fd02725b32eb55653f9ca25984dd98ad6",
    PREFIX + "true_Rm_patch_finite_N_C5_functions.py": "c4445ef8acfb33d29cf9022862f8b6396e6aa2f38bda7467409b5378a5deaf15",
    PREFIX + "true_Rm_patch_finite_N_C5_functions_check.json": "efd3b9a6fdfd4522ead1806d4ef4250864ef670909ca047d0fa85d9ce23bf207",
    PREFIX + "true_Rm_patch_C6_functions.py": "a6f1cb5deaf17e9482aa964aef7867aadbc815e4d703a415e4945cf665f14002",
    PREFIX + "true_Rm_patch_C6_functions_check.json": "0190fcf738fb8add5be85cf0767464a39cc0d18240e6b36ad4ec7930dce38aa7",
}

_OPTIONAL_PATCH_RECORDS = frozenset({
    PREFIX + "true_Rm_patch_finite_N_C5_functions.py",
    PREFIX + "true_Rm_patch_finite_N_C5_functions_check.json",
    PREFIX + "true_Rm_patch_C6_functions.py",
    PREFIX + "true_Rm_patch_C6_functions_check.json",
})

FIRST_NINE_CHARTS = (
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
POSTRESTORE_CHART = "postrestore"
PATCH_CHART = "actual_patch"
ALL_SOURCE_ROWS = FIRST_NINE_CHARTS + (POSTRESTORE_CHART, PATCH_CHART)

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
_ALL_N = PREFIX + "whole_Z_all_N_source_functions.py"
_ALL_N_CHECK = PREFIX + "whole_Z_all_N_source_functions_check.json"
_SWITCH = PREFIX + "whole_Z_all_N_switch_functions.py"
_SWITCH_CHECK = PREFIX + "whole_Z_all_N_switch_functions_check.json"
_LONG = PREFIX + "whole_Z_all_N_long_patch_functions.py"
_LONG_CHECK = PREFIX + "whole_Z_all_N_long_patch_functions_check.json"
_PRESCRIBED_BRIDGE_C6 = PREFIX + "prescribed_bridge_C6_functions.py"
_PRESCRIBED_BRIDGE_C6_CHECK = PREFIX + "prescribed_bridge_C6_functions_check.json"
_PRESCRIBED_R100_RM_C6 = PREFIX + "prescribed_R100_Rm_C6_functions.py"
_PRESCRIBED_R100_RM_C6_CHECK = PREFIX + "prescribed_R100_Rm_C6_functions_check.json"
_PATCH_C5 = PREFIX + "true_Rm_patch_finite_N_C5_functions.py"
_PATCH_C5_CHECK = PREFIX + "true_Rm_patch_finite_N_C5_functions_check.json"
_PATCH_C6 = PREFIX + "true_Rm_patch_C6_functions.py"
_PATCH_C6_CHECK = PREFIX + "true_Rm_patch_C6_functions_check.json"


def canonical_digest(value: Any) -> str:
    """Hash a graph without interpreting or evaluating any node operation."""

    raw = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _module_hash() -> str | None:
    try:
        return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    except OSError:
        return None


def _gate(gates: list[dict[str, Any]], name: str, **details: Any) -> None:
    row: dict[str, Any] = {"gate": name}
    row.update(details)
    if row not in gates:
        gates.append(row)


def _record(
    source_records: Mapping[str, Any],
    name: str,
    gates: list[dict[str, Any]],
    *,
    optional: bool = False,
) -> str | None:
    """Authenticate one text record; never read a path from the filesystem."""

    record = source_records.get(name) if isinstance(source_records, Mapping) else None
    if not isinstance(record, Mapping):
        if not optional:
            _gate(gates, "source_record_missing", record=name)
        return None
    expected_path = "experiments/root_st073/" + name
    text = record.get("source_text")
    actual = hashlib.sha256(text.encode("utf-8")).hexdigest() if isinstance(text, str) else None
    expected = RECORD_HASHES[name]
    if record.get("path") != expected_path:
        _gate(gates, "source_record_path_changed", record=name, expected=expected_path)
        return None
    if record.get("sha256") != expected or actual != expected:
        _gate(gates, "source_record_hash_changed", record=name, expected_sha256=expected)
        return None
    return text


def _parse(text: str, name: str, gates: list[dict[str, Any]]) -> ast.AST | None:
    try:
        return ast.parse(text, filename=name)
    except SyntaxError as exc:
        _gate(gates, "source_ast_invalid", record=name, error=str(exc))
        return None


def _definitions(
    tree: ast.AST | None,
    name: str,
    wanted: tuple[str, ...],
    gates: list[dict[str, Any]],
) -> dict[str, ast.AST]:
    found: dict[str, ast.AST] = {}
    if tree is not None:
        for item in ast.walk(tree):
            if isinstance(item, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)) and item.name in wanted:
                found.setdefault(item.name, item)
    for declaration in wanted:
        if declaration not in found:
            _gate(gates, "source_declaration_missing", record=name, declaration=declaration)
    return found


def _tokens(text: str, name: str, required: tuple[str, ...], gates: list[dict[str, Any]]) -> bool:
    missing = [token for token in required if token not in text]
    if missing:
        _gate(gates, "source_constructor_token_missing", record=name, tokens=missing)
        return False
    return True


def _method_text(tree: ast.AST | None, name: str, wanted: tuple[str, ...], gates: list[dict[str, Any]]) -> str:
    defs = _definitions(tree, name, wanted, gates)
    if not defs:
        return ""
    return "\n".join(ast.unparse(value) for value in defs.values())


def _name(node: ast.AST | None) -> str | None:
    return node.id if isinstance(node, ast.Name) else None


def _dict_field(function: ast.AST, variable: str, field: str) -> ast.AST | None:
    for node in ast.walk(function):
        if not isinstance(node, ast.Assign) or len(node.targets) != 1 or _name(node.targets[0]) != variable:
            continue
        value = node.value
        if isinstance(value, ast.Dict):
            for key, result in zip(value.keys, value.values):
                if isinstance(key, ast.Constant) and key.value == field:
                    return result
        elif isinstance(value, ast.Call) and isinstance(value.func, ast.Name) and value.func.id == "dict":
            for keyword in value.keywords:
                if keyword.arg == field:
                    return keyword.value
    return None


def _subscript_string(node: ast.AST | None) -> tuple[str | None, str | None]:
    if not isinstance(node, ast.Subscript) or not isinstance(node.value, ast.Name):
        return None, None
    slice_node = node.slice
    if isinstance(slice_node, ast.Constant) and isinstance(slice_node.value, str):
        return node.value.id, slice_node.value
    return node.value.id, None


def _c2_formula_witness(tree: ast.AST | None, name: str, gates: list[dict[str, Any]]) -> dict[str, Any]:
    functions = _definitions(tree, name, ("coefficient_C2", "build"), gates)
    coefficient = functions.get("coefficient_C2")
    build = functions.get("build")
    first_m = _name(_dict_field(coefficient, "first", "m")) if coefficient is not None else None
    second_m = _name(_dict_field(coefficient, "second", "m")) if coefficient is not None else None
    density_text = ast.unparse(build) if build is not None else ""
    density_tokens = (
        "coefficients[-1]",
        "coefficients[-2]",
        "invN",
        "alg.add",
        "alg.mul",
    )
    density_ok = _tokens(density_text, name + ":build", density_tokens, gates)
    ok = first_m == "B" and second_m == "zero" and density_ok
    if first_m != "B":
        _gate(gates, "c2_first_m_is_not_same_B", record=name, observed=first_m)
    if second_m != "zero":
        _gate(gates, "c2_second_m_is_not_zero", record=name, observed=second_m)
    return {
        "passed": ok,
        "first_m": first_m,
        "second_m": second_m,
        "density_formula": "coefficients[-1][m] + invN * coefficients[-2][m]",
        "ordinary_C2_m": "B",
        "density_formula_tokens": density_tokens,
    }


def _c5_formula_witness(tree: ast.AST | None, name: str, gates: list[dict[str, Any]]) -> dict[str, Any]:
    functions = _definitions(tree, name, ("loop_functions", "generic_roots", "finite_N_functions"), gates)
    loop = functions.get("loop_functions")
    roots = functions.get("generic_roots")
    finite = functions.get("finite_N_functions")
    loop_text = ast.unparse(loop) if loop is not None else ""
    roots_text = ast.unparse(roots) if roots is not None else ""
    finite_text = ast.unparse(finite) if finite is not None else ""
    loop_ok = _tokens(
        loop_text,
        name + ":loop_functions",
        ("roots['a']", "roots['b']", "T1hat", "B =", "B_over_Pstar", "roots['E']", "'1/2'"),
        gates,
    )
    roots_ok = _tokens(
        roots_text,
        name + ":generic_roots",
        ("C =", "E_y", "B =", "V_y", "aa", "bb", "t0", "generic['p2']"),
        gates,
    )
    finite_ok = _tokens(
        finite_text,
        name + ":finite_N_functions",
        ("B_over_Pstar", "invN", "s =", "inc", "m=s"),
        gates,
    )
    m_value = _dict_field(finite, "inc", "m") if finite is not None else None
    s_value = None
    if finite is not None:
        for node in ast.walk(finite):
            if isinstance(node, ast.Assign) and len(node.targets) == 1 and _name(node.targets[0]) == "s":
                s_value = node.value
                break
    loop_scale = ""
    if loop is not None:
        for node in ast.walk(loop):
            if isinstance(node, ast.Assign) and len(node.targets) == 1 and _name(node.targets[0]) == "B":
                loop_scale = ast.unparse(node.value)
                break
    assignments = {target.id: node.value for node in ast.walk(finite) if isinstance(node, ast.Assign)
        for target in node.targets if isinstance(target, ast.Name)} if finite is not None else {}
    exact = lambda value, expression: value is not None and ast.dump(value, include_attributes=False) == ast.dump(
        ast.parse(expression, mode="eval").body, include_attributes=False)
    scale_ok = exact(s_value, "a.scale(loop['B_over_Pstar'],invN)") and exact(assignments.get("invN"),
        "g.quotient(g.one,N,'same exact positive integer N on every source chart')")
    loop_B_ok = loop_scale and ast.dump(ast.parse(loop_scale, mode="eval").body, include_attributes=False) == ast.dump(
        ast.parse("lazy(a.scale(a.mul(roots['E'],M),'1/2'))", mode="eval").body, include_attributes=False)
    if not scale_ok:
        _gate(gates, "c5_m_scale_does_not_use_B_over_Pstar_over_N", record=name)
    if _name(m_value) != "s":
        _gate(gates, "c5_inc_m_is_not_scaled_s", record=name, observed=_name(m_value))
    return {
        "passed": bool(loop_ok and roots_ok and finite_ok and scale_ok and loop_B_ok and _name(m_value) == "s"),
        "exact_N_division_AST_bound": scale_ok,
        "exact_lazy_B_product_AST_bound": bool(loop_B_ok),
        "loop_B_definition": loop_scale,
        "finite_N_m_definition": "s = B_over_Pstar / N; inc['m'] = s",
        "generic_root_paths": {
            "C": "E - 2 * E_y",
            "B": "2 * V_y",
            "a": "C / E",
            "b": "B / E",
            "t0": "-B / C",
            "p2": "generic['p2']",
        },
    }


def _receipt(
    text: str | None,
    name: str,
    required: tuple[tuple[str, ...], ...],
    gates: list[dict[str, Any]],
) -> bool:
    if text is None:
        return False
    try:
        value = json.loads(text)
    except json.JSONDecodeError as exc:
        _gate(gates, "accepted_source_receipt_invalid", record=name, error=str(exc))
        return False
    ok = value.get("all_passed") is True
    if not ok:
        _gate(gates, "accepted_source_receipt_header_missing", record=name)
    for path in required:
        current: Any = value
        for key in path:
            current = current.get(key) if isinstance(current, Mapping) else None
        if current is not True:
            _gate(gates, "accepted_source_receipt_flag_missing", record=name, path=".".join(path))
            ok = False
    family = value.get("source_family")
    if isinstance(family, Mapping):
        for key, expected in SOURCE_FAMILY.items():
            if family.get(key) != expected:
                _gate(gates, "source_family_changed", record=name, field=key, expected=expected)
                ok = False
    else:
        for key, expected in SOURCE_FAMILY.items():
            if key in value and value.get(key) != expected:
                _gate(gates, "source_family_changed", record=name, field=key, expected=expected)
                ok = False
    return ok


def _provider_witnesses(texts: Mapping[str, str], gates: list[dict[str, Any]]) -> dict[str, Any]:
    """Check source-provider ownership without constructing a provider."""

    correlated_text = texts.get(_CORRELATED, "")
    switch_text = texts.get(_SWITCH, "")
    long_text = texts.get(_LONG, "")
    all_n_text = texts.get(_ALL_N, "")
    correlated_ok = _tokens(
        correlated_text,
        _CORRELATED,
        (
            "class WholeZCorrelatedBridgeFunctions",
            "def source_query",
            "packet['original_generic_source']",
            "original_C_B_E_C0_Z",
            "self.identity",
        ),
        gates,
    )
    switch_ok = _tokens(
        switch_text,
        _SWITCH,
        (
            "class WholeZAllNSwitchFunctions",
            "def leading_packet",
            "source_identity=self.identity",
            "exact_common_P0_axial5=live.P0",
            "original_generic_source=recovered",
            "original_roots=roots['roots']",
            "actual_phase_Z_exact_zero=True",
            "identical_source_P0_tuples_rebound_to_live_descendant_object=True",
        ),
        gates,
    )
    long_ok = _tokens(
        long_text,
        _LONG,
        (
            "class WholeZAllNLongPatchFunctions",
            "def leading_packet",
            "op.flow is not live.flow",
            "op.c is not live.c",
            "equivalent_rows",
            "exact_common_P0_axial5=live.P0",
            "original_generic_source=recovered",
            "actual_phase_Z_exact_zero=True",
        ),
        gates,
    )
    all_n_ok = _tokens(
        all_n_text,
        _ALL_N,
        (
            "class WholeZAllNBridgeSourceFunctions",
            "def primitive_cell",
            "exact_common_P0_axial5",
            "actual_phase_Z_exact_zero",
            "source_identity",
            "N_scaled_density_C0_Z",
            "original_roots",
        ),
        gates,
    )
    return {
        "first_micro": {"owner": "WholeZCorrelatedBridgeFunctions", "entry": "source_query", "passed": correlated_ok},
        "second_micro": {"owner": "WholeZCorrelatedBridgeFunctions", "entry": "source_query", "passed": correlated_ok},
        "frozen_macro": {"owner": "WholeZCorrelatedBridgeFunctions", "entry": "source_query", "passed": correlated_ok},
        "first_switch": {"owner": "WholeZAllNSwitchFunctions", "entry": "leading_packet", "passed": switch_ok},
        "second_switch": {"owner": "WholeZAllNSwitchFunctions", "entry": "leading_packet", "passed": switch_ok},
        "post_power": {"owner": "WholeZAllNSwitchFunctions", "entry": "leading_packet", "passed": switch_ok},
        "long_reshape": {"owner": "WholeZAllNLongPatchFunctions", "entry": "leading_packet", "passed": long_ok},
        "reference": {"owner": "WholeZAllNLongPatchFunctions", "entry": "leading_packet", "passed": long_ok},
        "restoration": {"owner": "WholeZAllNLongPatchFunctions", "entry": "leading_packet", "passed": long_ok},
        "postrestore": {"owner": "WholeZAllNLongPatchFunctions", "entry": "leading_packet", "passed": long_ok},
        "actual_patch": {"owner": "WholeZAllNLongPatchFunctions", "entry": "leading_packet", "passed": long_ok},
        "all_N_source_owner": all_n_ok,
    }


def _c5_owner_transfer_witness(texts: Mapping[str, str], gates: list[dict[str, Any]]) -> dict[str, Any]:
    """Check the prescribed C5 owner chain which supplies the actual fields.

    The C5 loop source receives its chart fields from the prescribed
    R100--Rm owner.  That owner receives the same prescribed bridge C6 owner,
    and the C5 constructor records the chart-by-chart source objects.  These
    source bindings are separate from the generic ``E``/``V`` names and are
    required before a C2 inherited ``B`` can be compared with C5's
    ``B_over_Pstar``.
    """

    c5 = texts.get(_C5, "")
    bridge = texts.get(_PRESCRIBED_BRIDGE_C6, "")
    r100 = texts.get(_PRESCRIBED_R100_RM_C6, "")
    c5_ok = _tokens(
        c5,
        _C5 + ":owner_transfer",
        (
            "import lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions as previous",
            "type(background) is not previous.CurrentOriginalPrescribedR100RmC6Functions",
            "self.background=background",
            "self.Pstar=background.Pstar;self.P0=background.P0",
            "self.source_charts={",
            "'first_micro':bridge.micro",
            "'second_micro':bridge.micro",
            "'frozen_macro':bridge.macro",
            "'first_switch':background.charts['first_switch']",
            "'long_reshape':background.charts['long_reshape']",
            "'reference':background.charts['reference']",
            "'restoration':background.charts['restoration']",
            "'postrestore':background.charts['post_restore']",
            "roots=generic_roots(g,chart,self.Pstar)",
        ),
        gates,
    )
    r100_ok = _tokens(
        r100,
        _PRESCRIBED_R100_RM_C6 + ":owner_transfer",
        (
            "import lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions as previous",
            "type(bridge_owner) is not previous.CurrentOriginalPrescribedBridgeC6Functions",
            "self.bridge=bridge_owner",
            "self.P0=bridge_owner.P0",
            "self.charts={}",
            "all_six_nonzero_incoming_and_every_chart_memory_retained=True",
            "original_independent_P0_and_physical_Pstar_F0_factors_once=True",
        ),
        gates,
    )
    bridge_ok = _tokens(
        bridge,
        _PRESCRIBED_BRIDGE_C6 + ":owner_transfer",
        (
            "class CurrentOriginalPrescribedBridgeC6Functions",
            "self.P0=self.pressure7()",
            "generic_fields",
            "fields=dict(E=",
            "m=a.scale(J['M'],inverse)",
            "P0=a.trim(P0)",
            "true_prescribed_background_before_finite_N_modification=True",
        ),
        gates,
    )
    return {
        "passed": c5_ok and r100_ok and bridge_ok,
        "c5_owner": "CurrentOriginalPrescribedLoopC5Functions.background",
        "r100_owner": "CurrentOriginalPrescribedR100RmC6Functions.bridge",
        "bridge_owner": "CurrentOriginalPrescribedBridgeC6Functions",
        "chart_source_transfer": {
            "first_micro": "bridge.micro",
            "second_micro": "bridge.micro",
            "frozen_macro": "bridge.macro",
            "first_switch": "background.charts['first_switch']",
            "second_switch": "background.charts['second_switch']",
            "post_power": "background.charts['post_switch']",
            "long_reshape": "background.charts['long_reshape']",
            "reference": "background.charts['reference']",
            "restoration": "background.charts['restoration']",
            "postrestore": "background.charts['post_restore']",
        },
        "field_transfer": "C5 generic_roots reads each transferred chart fields/radial/generic source object",
    }


def _optional_patch_witness(texts: Mapping[str, str], gates: list[dict[str, Any]]) -> bool:
    c5_text = texts.get(_PATCH_C5)
    c6_text = texts.get(_PATCH_C6)
    if c5_text is None or c6_text is None:
        return False
    c5_ok = _tokens(
        c5_text,
        _PATCH_C5,
        (
            "class CurrentOriginalTrueRmPatchFiniteNC5Functions",
            "engine.loop_functions",
            "engine.finite_N_functions",
            "transport_log",
            "self.N,self.eta,self.dstar",
            "actual_same_N_Section11_loop_inverse_density_and_history_C5=True",
        ),
        gates,
    )
    c6_ok = _tokens(
        c6_text,
        _PATCH_C6,
        (
            "class CurrentOriginalTrueRmPatchC6Functions",
            "self.background=modified_owner.background",
            "self.Pstar=self.background.Pstar;self.P0=self.background.P0",
            "exact_source_correspondence=dict(",
            "normalized_raw_histories=",
            "original_current_profile_recipes_bound_by_checked_C6_parent=True",
        ),
        gates,
    )
    return c5_ok and c6_ok


def _source_witnesses(source_records: Mapping[str, Any], gates: list[dict[str, Any]]) -> dict[str, Any]:
    texts: dict[str, str] = {}
    for name in RECORD_HASHES:
        text = _record(source_records, name, gates, optional=name in _OPTIONAL_PATCH_RECORDS)
        if text is not None:
            texts[name] = text

    required_names = set(RECORD_HASHES) - set(_OPTIONAL_PATCH_RECORDS)
    missing_required = sorted(required_names - set(texts))
    if missing_required:
        return {
            "records_present": sorted(texts),
            "required_records_missing": missing_required,
            "c2": {"passed": False},
            "c5": {"passed": False},
            "providers": {},
            "optional_patch_source": False,
            "same_B_provider_binding": False,
        }

    trees = {name: _parse(text, name, gates) for name, text in texts.items() if name.endswith(".py")}
    c2 = _c2_formula_witness(trees.get(_C2_DENSITY), _C2_DENSITY, gates)
    c5 = _c5_formula_witness(trees.get(_C5), _C5, gates)

    phase_defs = _definitions(trees.get(_C2_PHASE), _C2_PHASE, ("primitive_C2", "source_packet", "leading_source"), gates)
    phase_text = "\n".join(ast.unparse(node) for node in phase_defs.values())
    phase_ok = _tokens(
        phase_text,
        _C2_PHASE,
        (
            "common_velocity_E_axial5",
            "common_velocity_V_axial5",
            "actual_nonzero_b_axial5",
            "actual_Delta_axial5",
            "full_signed_p2_axial4",
            "actual_a_axial5",
            "actual_t0_axial5",
            "exact_b_and_t0_zero",
            "BZZ =",
            "old['B']",
            "exact_common_P0_axial5",
            "source_identity",
        ),
        gates,
    )
    # This cross-source check is deliberately textual and structural.  It
    # binds the C2 inherited B jet to the C5 E*M/2 source recipe while the
    # provider witnesses bind the actual roots and P0 object; no label alias
    # or numerical equality is accepted.
    owner_transfer = _c5_owner_transfer_witness(texts, gates)
    conditional_formula_ok = phase_ok and c2["passed"] and c5["passed"]
    # Both owner chains are visible, but no per-quantity functional identity
    # connects the C2 live providers to the C5 prescribed C6 chart objects.
    # Their existence, common family and formula names do not supply it.
    same_B = False
    _gate(gates, "actual_C2_provider_C5_C6_field_function_identity_unbound",
        required="per-chart E, a, b, p2, t0, Delta, exact parameters and phase function identities")

    receipt_requirements = {
        _C2_PHASE_CHECK: (
            ("current_original_17_chart_source_Z2_rows_installed",),
            ("current_original_same_inverse_C2_phase_primitive_functions_installed",),
            ("upstream_six_chart_C2_provider_installed",),
        ),
        _C2_DENSITY_CHECK: (
            ("current_original_17_chart_C2_signed_density_functions_installed",),
            ("current_original_17_chart_C2_own_rate_transport_and_five_target_functions_installed",),
            ("same_actual_terminal_amplitude_P0_source_and_global_phase_retained",),
        ),
        _C5_CHECK: (
            ("checks", "actual_source", "passed"),
            ("checks", "actual_source", "genuine_C6_inputs_and_C5_generic_roots"),
            ("checks", "actual_source", "same_original_N_eta_dstar_and_global_phase"),
            ("checks", "actual_source", "all_ten_correction_predecessors_preserved"),
            ("current_original_true_prescribed_ten_chart_signed_loop_C5_functions_installed",),
            ("current_original_true_prescribed_ten_chart_finite_N_density_C5_functions_installed",),
        ),
        _BRIDGE_CHECK: (("current_original_whole_Z_live_leading_bridge_source_functions_installed",),),
        _CORRELATED_CHECK: (("current_original_whole_Z_correlated_bridge_C1_functions_and_Rc_retransport_installed",),),
        _ALL_N_CHECK: (("current_original_whole_Z_bridge_full_phase_all_N_Nscaled_C1_functions_and_R100_transport_installed",),),
        _SWITCH_CHECK: (("current_original_whole_Z_full_phase_all_N_switch_C1_functions_and_normalized_R110_transport_installed",),),
        _LONG_CHECK: (("current_original_whole_Z_full_phase_all_N_long_patch_C1_functions_and_normalized_Rh_transport_installed",),),
        _PRESCRIBED_BRIDGE_C6_CHECK: (
            ("current_original_true_comparison_bridge_C7_functions_installed",),
            ("current_original_true_prescribed_micro_macro_C6_functions_installed",),
            ("current_original_true_prescribed_bridge_Q_inertial_p2_C5_functions_installed",),
        ),
        _PRESCRIBED_R100_RM_C6_CHECK: (
            ("current_original_true_prescribed_R100_R110_C6_functions_installed",),
            ("current_original_true_prescribed_finite_long_reference_restore_Rm_C6_functions_installed",),
            ("current_original_true_prescribed_R100_Rm_Q_inertial_p2_C5_functions_installed",),
        ),
    }
    receipt_ok = {
        name: _receipt(texts.get(name), name, required, gates)
        for name, required in receipt_requirements.items()
    }
    providers = _provider_witnesses(texts, gates)
    for chart in ALL_SOURCE_ROWS:
        owner = providers.get(chart, {})
        owner_ok = isinstance(owner, Mapping) and owner.get("passed") is True
        if chart in FIRST_NINE_CHARTS and not owner_ok:
            _gate(gates, "source_provider_binding_unproved", chart=chart)
    provider_receipts_ok = all(receipt_ok.values())
    if not owner_transfer["passed"]:
        _gate(gates, "c5_prescribed_owner_transfer_unproved")
    optional_patch = _optional_patch_witness(texts, gates)
    if texts.get(_PATCH_C5_CHECK) is not None:
        patch_receipts_ok = _receipt(
            texts.get(_PATCH_C5_CHECK),
            _PATCH_C5_CHECK,
            (
                ("current_original_true_Rm_patch_finite_N_loop_density_history_C5_functions_installed",),
                ("current_original_true_Rm_patch_finite_N_C5_transport_installed",),
            ),
            gates,
        )
        optional_patch = optional_patch and patch_receipts_ok
    if texts.get(_PATCH_C6_CHECK) is not None:
        optional_patch = optional_patch and _receipt(
            texts.get(_PATCH_C6_CHECK),
            _PATCH_C6_CHECK,
            (
                ("current_original_true_Rm_leading_implicit_control_C6_functions_installed",),
                ("current_original_true_Rm_patch_background_C6_generic_C5_functions_installed",),
            ),
            gates,
        )
    return {
        "records_present": sorted(texts),
        "c2": c2,
        "c5": c5,
        "phase_source_roots": {
            "E": "common_velocity_E_axial5",
            "V": "common_velocity_V_axial5",
            "a": "actual_a_axial5 / actual_correlated_a_axial5 / full_original_shear_a_axial5",
            "b": "actual_nonzero_b_axial5, or exact_b_and_t0_zero",
            "t0": "actual_t0_axial5, or exact_b_and_t0_zero",
            "p2": "full_signed_p2_axial4 / exact_full_p2_axial4",
            "Delta": "actual_Delta_axial5 / Delta_axial5",
        },
        "providers": providers,
        "receipts": receipt_ok,
        "c5_prescribed_owner_transfer": owner_transfer,
        "conditional_density_linear_identity": bool(conditional_formula_ok),
        "same_B_provider_binding": same_B and provider_receipts_ok and all(
            isinstance(providers.get(chart), Mapping) and providers[chart].get("passed") is True
            for chart in FIRST_NINE_CHARTS
        ),
        "optional_patch_source": optional_patch,
    }


def _row(chart: str, admitted: bool, *, reason: str | None = None, provider: Mapping[str, Any] | None = None) -> dict[str, Any]:
    row = {
        "chart": chart,
        "status": "proved" if admitted else "unresolved",
        "exact_identity": "C2_m = B = N * (B_over_Pstar / N) = N * C5_m" if admitted else None,
        "c2_density_m": "coefficient_C2: first['m'] = B; second['m'] = zero",
        "c5_density_m": "finite_N_functions: s = B_over_Pstar / N; inc['m'] = s",
        "same_B_is_provider_bound": admitted,
        "provider": dict(provider) if isinstance(provider, Mapping) else None,
        "source_only": True,
        "numeric_values_computed": False,
        "operator_calls": 0,
        "jet_indices_promoted": False,
    }
    if reason is not None:
        row["unresolved_gate"] = reason
    return row


def derive(nodes: Sequence[Any], source_records: Mapping[str, Any], *, expected_nodes_digest: str) -> dict[str, Any]:
    """Return a fail-closed source certificate for the C2/C5 m rows."""

    if not isinstance(expected_nodes_digest, str) or expected_nodes_digest != EXPECTED_NODES_DIGEST:
        raise ValueError("unexpected source graph digest assertion")
    try:
        node_count = len(nodes)
        actual_digest = canonical_digest(nodes)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("nodes are not a canonical JSON graph") from exc
    if node_count != EXPECTED_NODE_COUNT or actual_digest != EXPECTED_NODES_DIGEST:
        raise ValueError("nodes canonical digest or node count does not match expected source graph")

    gates: list[dict[str, Any]] = []
    witness = _source_witnesses(source_records, gates)
    bridge_ok = witness.get("same_B_provider_binding") is True
    rows: list[dict[str, Any]] = []
    for chart in FIRST_NINE_CHARTS:
        provider = witness.get("providers", {}).get(chart)
        rows.append(
            _row(
                chart,
                bridge_ok and isinstance(provider, Mapping) and provider.get("passed") is True,
                reason=None if bridge_ok and isinstance(provider, Mapping) and provider.get("passed") is True else "first_n_same_B_source_provider_binding",
                provider=provider,
            )
        )
    post_provider = witness.get("providers", {}).get(POSTRESTORE_CHART)
    post_ok = bridge_ok and isinstance(post_provider, Mapping) and post_provider.get("passed") is True
    rows.append(_row(POSTRESTORE_CHART, post_ok, reason=None if post_ok else "postrestore_direct_source_provider", provider=post_provider))
    patch_provider = witness.get("providers", {}).get(PATCH_CHART)
    patch_ok = bridge_ok and witness.get("optional_patch_source") is True and isinstance(patch_provider, Mapping) and patch_provider.get("passed") is True
    rows.append(_row(PATCH_CHART, patch_ok, reason=None if patch_ok else "patch_direct_C5_source_theorem", provider=patch_provider))

    if not post_ok:
        _gate(gates, "postrestore_direct_source_provider_unproved")
    if not patch_ok:
        _gate(gates, "patch_direct_C5_source_theorem_unproved")
    if not bridge_ok:
        _gate(gates, "same_B_provider_chain_unproved")

    admitted_rows = [row for row in rows if row["status"] == "proved"]
    unresolved_rows = [row for row in rows if row["status"] != "proved"]
    return {
        "schema": SCHEMA,
        "source_only": True,
        "passed": bool(admitted_rows) and not gates,
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
        "rows": rows,
        "admitted_rows": admitted_rows,
        "unresolved_rows": unresolved_rows,
        "unresolved_gates": gates,
        "component_proofs": {
            "c2_formula": witness.get("c2", {}),
            "c5_formula": witness.get("c5", {}),
            "same_B_provider_binding": witness.get("same_B_provider_binding", False),
            "conditional_density_linear_identity": witness.get("conditional_density_linear_identity", False),
            "c5_prescribed_owner_transfer": witness.get("c5_prescribed_owner_transfer", {}),
            "source_provider_bindings": witness.get("providers", {}),
            "source_receipts": witness.get("receipts", {}),
            "source_records_present": witness.get("records_present", []),
        },
        "source_root_paths": witness.get("phase_source_roots", {}),
        "preserved_source_handle_assertions": {
            "N_handle": 18,
            "eta_handle": 21,
            "d_star_handle": 22,
            "handle_numbers_are_not_parameter_values": True,
            "cross_owner_constant_equality_proved": False,
            "phase_origin": "same original global phase origin; no new patch phase origin",
        },
        "numeric_values_computed": False,
        "mathematical_operator_replayed": False,
        "old_operators_replayed": False,
        "operator_calls": 0,
        "producer_provider_builder_checker_replayed": False,
        "runtime_metadata": None,
        "record_hashes": dict(RECORD_HASHES),
        "module_sha256": _module_hash(),
        "proof_limitations": [
            "Actual C2/C5 source function equality remains unproved; all actual density rows remain unadmitted.",
            "The conditional algebra is C2_m=B and C5_m=B_over_Pstar/N; it is not a proof that their B functions agree.",
            "The certificate checks source text, AST structure, and accepted receipts only.",
            "No saved graph node is interpreted as a mathematical operator.",
            "A missing direct patch C5 source receipt leaves actual_patch unresolved.",
            "No quantitative bounds, jet indices, guard strings, or numerical root values are promoted.",
        ],
    }


__all__ = ["RECORD_HASHES", "canonical_digest", "derive"]
