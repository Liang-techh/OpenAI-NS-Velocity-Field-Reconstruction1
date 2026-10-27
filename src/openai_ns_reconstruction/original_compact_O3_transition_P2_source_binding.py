"""Source-binding-only P2 kernel draft (v3).

This file is a static binder.  It accepts the frozen manager wrapper, unwraps
``source['input']`` once, hashes the supplied full source text, parses ASTs,
and checks bounded literal records.  It never imports or calls a project
provider, checker, normalizer, reducer, quadrature, recurrence, field,
integral, or proof implementation.
"""

from __future__ import annotations

import ast
import hashlib
from typing import Any, Callable, Mapping


INPUT_SCHEMA = "compact_O3_transition_P2_source_binding_input_v1"
RESULT_SCHEMA = "original_compact_O3_transition_P2_source_binding"
COMPONENT = "compact_O3_transition_P2_source_binding"
TRACE_SHA = "b625a9340e914f3ca91e3418288338cb70d868aa54894e3984617679c51a888e"
NATIVE_ARCHIVE_SHA = "75895e1637411cfd2e490778f495f3adf904295a8840184cdc467e006b23305d"
C3_ARCHIVE_SHA = "cf51e2131c3d1cb76f93a2deb76af12421e2deb75f7888bd310084b3b3eec193"

NATIVE_ROOT_IDS = (242304, 242308, 242314, 242321)
C3_ROOT_IDS = (4598, 4599, 9992, 21729)
ORDINARY_Z = (0, 1, 2, 3)
FACTORIALS = {0: 1, 1: 1, 2: 2, 3: 6}

NATIVE_FORMULA = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py"
NATIVE_ROOTS = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
NATIVE_BACKGROUND = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.py"
C3_RECOVERY = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rm_generic_inputs.py"
C3_REFERENCE = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_reference_restore_finite_N.py"
C3_LONG = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_long_reshape_finite_N.py"
C3_O3 = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_O3_Rc_finite_N.py"
C3_WHOLE = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_outer_Rc_functions.py"
C3_WHOLE_O3 = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_whole_Z_O3_Rc_finite_N.py"
C3_REPAIR = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_C3_repair_band.py"
C2_PROJECTION = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_C2_phase_primitives.py"
C3_PROJECTION = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_C3_source_targets.py"

SOURCE_SHA = {
    NATIVE_BACKGROUND: "4464a934d303a95a0da157e9e289f6ed8a04d5ff0adcb0d663857f0d9c57dcd3",
    NATIVE_FORMULA: "7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7",
    NATIVE_ROOTS: "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d",
    C3_RECOVERY: "62df6a6513b6e9c521cb6cf7cff25710dec17bf04f7b23f4bc0f15dbc1a57a88",
    C3_REFERENCE: "f338527b4139f549f2f76faaaa65f94923c423bfbc0d133151a7aca0dfa3a47d",
    C3_LONG: "07477c8d4b7e5b7e25993415ba55be0ddaf0ccff1b504045a55d5aa7bf3612c4",
    C3_O3: "45284dea25aa4c2ff06f2fb925e7efd85547d18b38c106d39419618555756a48",
    C3_WHOLE: "9d73659a1c76163dddee6653af53f56ebffd7045986d871dc559ba1e053620b7",
    C3_WHOLE_O3: "83e5a8bdae75d20a261d24dd22d30cc7cb5e603e8fa2f0133958ca880b5e76db",
    C3_REPAIR: "62e2fa8598bf307126be0ffe0b4e23de83152104524db03f61b4944c225132e4",
    C2_PROJECTION: "84f5da3eb64096f82df88dd828b22b09bce18192a03925d91cfee0b193902d2c",
    C3_PROJECTION: "de41858de44ff9cea9b5d2b17aeca55eeca2c114969ca4df383ec1345b6badeb",
}
REQUIRED_SOURCES = tuple(SOURCE_SHA)

NEW_FLAGS = (
    "O3_transition_P2_source_ports_bound",
    "O3_transition_P2_pressure_restoration_route_bound",
    "O3_transition_P2_ordinary_Z_source_projection_bound",
)
COUNT_KEYS = (
    "source_port_rows",
    "native_source_port_rows",
    "C3_source_port_rows",
    "projection_rows",
    "pressure_route_sides",
    "new_mathematical_proofs",
    "new_algebra_family_calls",
    "numeric_calls",
    "source_provider_calls",
    "old_math_calls",
)


class _Failure(Exception):
    pass


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _ast_sha(node: ast.AST) -> str:
    return _sha(ast.dump(node, annotate_fields=True, include_attributes=False).encode("utf-8"))


def _text(node: ast.AST) -> str:
    try:
        return ast.unparse(node)
    except Exception:
        return ast.dump(node, annotate_fields=True, include_attributes=False)


def _is_name(node: Any, value: str) -> bool:
    return isinstance(node, ast.Name) and node.id == value


def _is_const(node: Any, value: Any) -> bool:
    return isinstance(node, ast.Constant) and node.value == value


def _constants(node: ast.AST) -> list[Any]:
    return [item.value for item in ast.walk(node) if isinstance(item, ast.Constant)]


def _names(node: ast.AST) -> set[str]:
    return {item.id for item in ast.walk(node) if isinstance(item, ast.Name)}


def _attr_chain(node: Any) -> list[str] | None:
    parts: list[str] = []
    while isinstance(node, ast.Attribute):
        parts.append(node.attr)
        node = node.value
    if not isinstance(node, ast.Name):
        return None
    parts.append(node.id)
    return list(reversed(parts))


def _is_call_chain(node: Any, chain: list[str]) -> bool:
    return isinstance(node, ast.Call) and _attr_chain(node.func) == chain


def _subscript(node: Any, base: str, key: Any | None = None) -> bool:
    if not isinstance(node, ast.Subscript) or not _is_name(node.value, base):
        return False
    if key is None:
        return True
    return _is_const(node.slice, key)


def _find(tree: ast.AST, name: str, owner: str | None = None) -> ast.AST | None:
    for node in ast.walk(tree):
        if owner is None and isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name:
            return node
        if owner is not None and isinstance(node, ast.ClassDef) and node.name == owner:
            for child in node.body:
                if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)) and child.name == name:
                    return child
    return None


def _assignment_nodes(node: ast.AST, name: str) -> list[ast.Assign | ast.AnnAssign]:
    result: list[ast.Assign | ast.AnnAssign] = []
    for item in ast.walk(node):
        targets = item.targets if isinstance(item, ast.Assign) else [item.target] if isinstance(item, ast.AnnAssign) else []
        if any(_is_name(target, name) for target in targets):
            result.append(item)
    return result


def _subscript_assignment_nodes(node: ast.AST, base: str) -> list[ast.Assign | ast.AnnAssign]:
    result: list[ast.Assign | ast.AnnAssign] = []
    for item in ast.walk(node):
        targets = item.targets if isinstance(item, ast.Assign) else [item.target] if isinstance(item, ast.AnnAssign) else []
        if any(_subscript(target, base) for target in targets):
            result.append(item)
    return result


def _dict_keyword(call: Any, key: str) -> ast.keyword | None:
    if not isinstance(call, ast.Call):
        return None
    for item in call.keywords:
        if item.arg == key:
            return item
    return None


def _literal_nodes(node: ast.AST) -> list[dict[str, Any]]:
    return [
        {"lineno": getattr(item, "lineno", None), "value_repr": repr(item.value)}
        for item in ast.walk(node)
        if isinstance(item, ast.Constant)
        and (isinstance(item.value, (str, int, float, bool)) or item.value is None)
    ]


def _selector_ast_audit(tree: ast.AST, selectors: Any) -> list[dict[str, Any]]:
    """Bind selector metadata to actual AST nodes where the selector is typed.

    One historical trace selector has a truncated AST digest.  The complete
    source-record function digest remains authoritative; the malformed trace
    digest is retained as an audit item instead of becoming a text guard.
    """

    result = []
    if not isinstance(selectors, Mapping):
        return result
    for selector, meta in selectors.items():
        if not isinstance(meta, Mapping) or "ast_sha256" not in meta:
            continue
        owner = None
        name = selector
        if "." in selector:
            owner, name = selector.rsplit(".", 1)
        node = _find(tree, name, owner)
        expected = meta.get("ast_sha256")
        actual = _ast_sha(node) if node is not None else None
        result.append(
            {
                "selector": selector,
                "trace_ast_sha256": expected,
                "actual_ast_sha256": actual,
                "trace_digest_well_formed": isinstance(expected, str) and len(expected) == 64,
                "actual_selector_bound": node is not None,
                "matches": actual == expected,
            }
        )
    return result


def _bind_source(record: Any, path: str) -> dict[str, Any]:
    if not isinstance(record, Mapping) or record.get("path") != path:
        raise _Failure(f"fixed source path mismatch: {path}")
    text = record.get("source_text")
    if not isinstance(text, str):
        raise _Failure(f"full source text missing: {path}")
    raw = text.encode("utf-8")
    if record.get("bytes") != len(raw):
        raise _Failure(f"source byte count mismatch: {path}")
    actual = _sha(raw)
    if actual != SOURCE_SHA[path] or record.get("sha") != actual:
        raise _Failure(f"fixed source SHA mismatch: {path}")
    tree = ast.parse(text, filename=path)
    actual_ast = _ast_sha(tree)
    if record.get("module_ast_sha256") != actual_ast:
        raise _Failure(f"module AST SHA mismatch: {path}")
    function_records = record.get("functions")
    if not isinstance(function_records, Mapping):
        raise _Failure(f"source function metadata missing: {path}")
    selectors = record.get("selectors_from_trace")
    if not isinstance(selectors, Mapping):
        raise _Failure(f"source selector metadata missing: {path}")
    for name, info in function_records.items():
        node = _find(tree, name)
        if node is not None and isinstance(info, Mapping) and info.get("ast_sha256") != _ast_sha(node):
            raise _Failure(f"function AST SHA mismatch: {path}:{name}")
    return {
        "path": path,
        "sha": actual,
        "module_ast_sha256": actual_ast,
        "tree": tree,
        "selector_ast_audit": _selector_ast_audit(tree, selectors),
    }


def _emit(events: list[dict[str, Any]], journal: Any, event: dict[str, Any]) -> None:
    events.append(event)
    if journal is not None:
        if not callable(journal):
            raise _Failure("journal must be callable")
        journal(event)


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise _Failure(f"missing mapping: {label}")
    return value


def _root_records(slice_data: Mapping[str, Any], ids: tuple[int, ...], digest: str, side: str) -> dict[int, dict[str, Any]]:
    if slice_data.get("archive_sha256") != digest:
        raise _Failure(f"{side} archive SHA mismatch")
    if slice_data.get("graph_key") != "exact_graph_nodes" or slice_data.get("recursive_graph_walk") is not False:
        raise _Failure(f"{side} graph metadata changed")
    expected_length = 258857 if side == "native" else 24952
    if slice_data.get("decoded_graph_length") != expected_length:
        raise _Failure(f"{side} decoded graph length missing")
    if slice_data.get("selection_policy") != "four explicit roots plus one immediate direct graph level":
        raise _Failure(f"{side} graph selection policy changed")
    if tuple(slice_data.get("root_ids", ())) != ids:
        raise _Failure(f"{side} root IDs changed")
    roots = slice_data.get("roots")
    if not isinstance(roots, Mapping) or set(roots) != {str(item) for item in ids}:
        raise _Failure(f"{side} root map is not exact")
    result = {}
    for node_id in ids:
        row = roots[str(node_id)]
        if not isinstance(row, Mapping) or row.get("node") != node_id or not isinstance(row.get("record"), Mapping):
            raise _Failure(f"{side} root record malformed: {node_id}")
        result[node_id] = dict(row["record"])
    return result


def _native_literals(slice_data: Mapping[str, Any]) -> dict[str, Any]:
    roots = _root_records(slice_data, NATIVE_ROOT_IDS, NATIVE_ARCHIVE_SHA, "native")
    expected_roots = {
        242304: {"operation": "positive_quotient", "numerator": 242298, "denominator": 241552, "source_positive_certificate": "same actual positive modified swirl E_N"},
        242308: {"operation": "positive_quotient", "numerator": 242307, "denominator": 241552, "source_positive_certificate": "same actual positive modified swirl E_N"},
        242314: {"operation": "positive_quotient", "numerator": 242313, "denominator": 241552, "source_positive_certificate": "same actual positive modified swirl E_N"},
        242321: {"operation": "positive_quotient", "numerator": 242320, "denominator": 241552, "source_positive_certificate": "same actual positive modified swirl E_N"},
    }
    if roots != expected_roots:
        raise _Failure("native P2 root records changed")
    direct = slice_data.get("direct_records")
    if not isinstance(direct, Mapping) or set(direct) != {"241552", "242298", "242307", "242313", "242320"}:
        raise _Failure("native direct record set changed")
    expected_direct = {
        "241552": {"operation": "product", "arguments": [233688, 241551]},
        "242298": {"operation": "product", "arguments": [242292, 241787]},
        "242307": {"operation": "sum", "arguments": [242299, 242306]},
        "242313": {"operation": "sum", "arguments": [242300, 242312]},
        "242320": {"operation": "sum", "arguments": [242301, 242319]},
    }
    for key, expected in expected_direct.items():
        row = direct[key]
        if not isinstance(row, Mapping) or row.get("node") != int(key) or row.get("record") != expected:
            raise _Failure(f"native direct record changed: {key}")
    return {
        "side": "native",
        "archive_path": slice_data.get("archive_path"),
        "archive_sha256": slice_data.get("archive_sha256"),
        "chart": "O3_transition",
        "coordinate_node": 241536,
        "root_ids": list(NATIVE_ROOT_IDS),
        "p2_root_selectors": [
            {
                "order": order,
                "root_node": node_id,
                "quotient_root_selector": True,
                "root_numerator_is_raw_W": False,
                "p2_numerator_role": (
                    "Z0 quotient numerator; not raw W0"
                    if order == 0
                    else "Z-order recurrence-adjusted quotient numerator; not raw W_j"
                ),
                "dependency_order": order + 1,
            }
            for order, node_id in enumerate(NATIVE_ROOT_IDS)
        ],
        "positive_E_denominator_node": 241552,
        "bounded_direct_records": {key: direct[key]["record"] for key in sorted(direct)},
        "raw_W_graph_correspondence": {"W0_node": 242298, "is_authoritative_W_source": False},
    }


def _c3_literals(slice_data: Mapping[str, Any]) -> dict[str, Any]:
    roots = _root_records(slice_data, C3_ROOT_IDS, C3_ARCHIVE_SHA, "rebuilt C3")
    projection = []
    for order, node_id in enumerate(C3_ROOT_IDS):
        row = roots[node_id]
        recipe = row.get("recipe") if isinstance(row.get("recipe"), Mapping) else None
        paths = recipe.get("quantity_paths") if isinstance(recipe, Mapping) and isinstance(recipe.get("quantity_paths"), Mapping) else None
        binding = row.get("source_projection_binding")
        if (
            row.get("operation") != "current_original_leading_function_recipe"
            or row.get("quantity") != "p2"
            or row.get("native_chart") != "O3_transition"
            or row.get("coordinate") != 4582
            or row.get("Z_variable") != 20
            or row.get("Z_order") != order
            or recipe is None
            or recipe.get("module") != C3_WHOLE.rsplit("/", 1)[-1]
            or recipe.get("function") != "WholeZAllNOuterRcFunctions.leading_packet"
            or recipe.get("defining_provider") != "WholeZAllNOuterRcFunctions.leading_packet"
            or paths is None
            or paths.get("p2") != "original_roots.p2[(0,Z_order)]"
        ):
            raise _Failure(f"rebuilt C3 recipe metadata changed: {node_id}")
        if order >= 2 and row.get("derivative_of_range_endpoint") is not False:
            raise _Failure(f"rebuilt C3 derivative endpoint marker changed: {node_id}")
        if order < 2 and "derivative_of_range_endpoint" in row:
            raise _Failure(f"rebuilt C3 low-order endpoint marker changed: {node_id}")
        if order >= 2 and row.get("Taylor_coefficient_factorial") != FACTORIALS[order]:
            raise _Failure(f"rebuilt C3 factorial metadata changed: {node_id}")
        if order == 2:
            if not isinstance(binding, Mapping) or binding.get("module") != C2_PROJECTION.rsplit("/", 1)[-1] or binding.get("module_sha256") != SOURCE_SHA[C2_PROJECTION] or binding.get("function") != "CurrentC2PhasePrimitives.source_packet" or binding.get("AST_sha256") != "3433d95a19fcae8105f7654f5fd6618473c365f6592ff23c0b936813d2c48917":
                raise _Failure("C3 order-2 source projection binding changed")
        if order == 3:
            if not isinstance(binding, Mapping) or binding.get("module") != C3_PROJECTION.rsplit("/", 1)[-1] or binding.get("module_sha256") != SOURCE_SHA[C3_PROJECTION] or binding.get("function") != "CurrentC3SourceTargets.source_packet" or binding.get("AST_sha256") != "8a994109bbb92735130640bd3564ad23bcf95daa7f4d346ea86fb6905b681338":
                raise _Failure("C3 order-3 source projection binding changed")
        projection.append(
            {
                "order": order,
                "recipe_node": node_id,
                "quotient_root_selector": True,
                "root_numerator_is_raw_W": False,
                "dependency_order": order + 1,
                "Taylor_coefficient_factorial": FACTORIALS[order],
                "factorial_applied_once": order in (2, 3),
                "source_projection_binding": binding,
                "quantity_path": paths.get("p2"),
            }
        )
    return {
        "side": "rebuilt_C3",
        "archive_path": slice_data.get("archive_path"),
        "archive_sha256": slice_data.get("archive_sha256"),
        "chart": "O3_transition",
        "coordinate_node": 4582,
        "Z_variable_node": 20,
        "root_ids": list(C3_ROOT_IDS),
        "p2_root_selectors": projection,
        "factorial_map": {str(key): value for key, value in FACTORIALS.items()},
    }


def _native_formula(sources: Mapping[str, Any]) -> dict[str, Any]:
    source = sources[NATIVE_FORMULA]
    roots = sources[NATIVE_ROOTS]
    fn = _find(source["tree"], "generic_functions")
    roots_fn = _find(roots["tree"], "generic_roots")
    if fn is None or roots_fn is None:
        raise _Failure("native formula creator missing")
    iz_nodes = _assignment_nodes(fn, "Iz")
    p2_nodes = _assignment_nodes(fn, "p2")
    if len(iz_nodes) != 1 or len(p2_nodes) != 1:
        raise _Failure("native Iz/p2 assignment multiplicity changed")
    iz = iz_nodes[0].value
    p2 = p2_nodes[0].value
    if not (_is_call_chain(iz, ["a", "scale"]) and len(iz.args) == 2 and _is_name(iz.args[1], "radius")):
        raise _Failure("native Iz assignment is not scale(..., radius)")
    if not _is_call_chain(iz.args[0], ["a", "add"]) or len(iz.args[0].args) != 2:
        raise _Failure("native Iz sum shape changed")
    if not _is_name(iz.args[0].args[0], "izl") or not (
        _is_call_chain(iz.args[0].args[1], ["a", "scale"])
        and len(iz.args[0].args[1].args) == 2
        and _is_name(iz.args[0].args[1].args[0], "izq")
        and _is_name(iz.args[0].args[1].args[1], "Pstar")
    ):
        raise _Failure("native Iz linear/quadratic operand shape changed")
    if not (_is_call_chain(p2, ["a", "div"]) and len(p2.args) >= 2 and _is_name(p2.args[0], "Iz") and _is_name(p2.args[1], "E")):
        raise _Failure("native p2 quotient operand shape changed")
    root_p2 = _assignment_nodes(roots_fn, "p2")
    if len(root_p2) != 1 or not isinstance(root_p2[0].value, ast.Subscript) or not _is_name(root_p2[0].value.value, "generic") or not _is_const(root_p2[0].value.slice, "p2"):
        raise _Failure("native generic_roots p2 handoff shape changed")
    return {
        "source_path": source["path"],
        "source_sha": source["sha"],
        "function": "generic_functions",
        "function_ast_sha256": _ast_sha(fn),
        "assignment_target": "Iz",
        "assignment_ast_sha256": _ast_sha(iz_nodes[0]),
        "assignment_source": _text(iz_nodes[0]),
        "W_side": "R*(I_z_linear + Pstar*I_z_quadratic)",
        "P2_side": "W_side/E",
        "root_handoff": {
            "source_path": roots["path"],
            "source_sha": roots["sha"],
            "function": "generic_roots",
            "function_ast_sha256": _ast_sha(roots_fn),
            "assignment_ast_sha256": _ast_sha(root_p2[0]),
        },
    }


def _exact_histories_assignment(fn: ast.AST) -> ast.Assign:
    rows = _assignment_nodes(fn, "histories")
    if len(rows) != 1 or not isinstance(rows[0].value, ast.DictComp):
        raise _Failure("C3 histories assignment shape changed")
    comp = rows[0].value
    if not isinstance(comp.value, ast.IfExp) or not _is_call_chain(comp.value.body, ["scale"]):
        raise _Failure("C3 invS conditional shape changed")
    scaled = comp.value.body
    if len(scaled.args) != 2 or not (_subscript(scaled.args[0], "rows", 0) and _is_name(scaled.args[1], "invS")):
        raise _Failure("C3 invS branch operand shape changed")
    if not _subscript(comp.value.orelse, "rows", 0):
        raise _Failure("C3 unscaled history branch changed")
    if not isinstance(comp.value.test, ast.Compare) or len(comp.value.test.ops) != 1 or not isinstance(comp.value.test.ops[0], ast.In) or not _is_name(comp.value.test.left, "key"):
        raise _Failure("C3 m/k history predicate shape changed")
    comparator = comp.value.test.comparators[0]
    if not isinstance(comparator, ast.Tuple) or [_constant_value(x) for x in comparator.elts] != ["m", "k"]:
        raise _Failure("C3 m/k history predicate values changed")
    if not _is_name(comp.key, "key") or len(comp.generators) != 1 or not isinstance(comp.generators[0].target, ast.Tuple):
        raise _Failure("C3 history generator shape changed")
    target = comp.generators[0].target
    history_iter = comp.generators[0].iter
    if not (
        isinstance(history_iter, ast.Call)
        and isinstance(history_iter.func, ast.Attribute)
        and history_iter.func.attr == "items"
        and _subscript(history_iter.func.value, "raw", "histories")
    ):
        raise _Failure("C3 history source packet shape changed")
    if len(target.elts) != 2 or not (_is_name(target.elts[0], "key") and _is_name(target.elts[1], "rows")):
        raise _Failure("C3 history source packet shape changed")
    exact_scale_count = sum(
        isinstance(item, ast.Call)
        and _is_call_chain(item, ["scale"])
        and len(item.args) == 2
        and _subscript(item.args[0], "rows", 0)
        and _is_name(item.args[1], "invS")
        for item in ast.walk(fn)
    )
    if exact_scale_count != 1:
        raise _Failure("C3 m/k invS normalization is not applied exactly once")
    return rows[0]


def _constant_value(node: Any) -> Any:
    return node.value if isinstance(node, ast.Constant) else object()


def _c3_provider_return(fn: ast.AST) -> dict[str, Any]:
    result_nodes = _assignment_nodes(fn, "result")
    if len(result_nodes) != 1 or not _is_name(result_nodes[0].value.func if isinstance(result_nodes[0].value, ast.Call) else None, "dict"):
        raise _Failure("C3 result dictionary assignment missing")
    result_call = result_nodes[0].value
    generic_kw = _dict_keyword(result_call, "actual_generic_source_numerators")
    sectors_kw = _dict_keyword(result_call, "full_signed_inertial_sectors_axial4")
    if generic_kw is None or sectors_kw is None:
        raise _Failure("C3 raw W return keys missing")
    if not isinstance(generic_kw.value, ast.Call) or not _is_name(generic_kw.value.func, "dict") or not isinstance(sectors_kw.value, ast.Call) or not _is_name(sectors_kw.value.func, "dict"):
        raise _Failure("C3 raw W return dictionaries changed")
    inertial_kw = _dict_keyword(generic_kw.value, "inertial_axial")
    linear_kw = _dict_keyword(sectors_kw.value, "axial_linear")
    quadratic_kw = _dict_keyword(sectors_kw.value, "axial_quadratic")
    if inertial_kw is None or not _is_name(inertial_kw.value, "Iz") or linear_kw is None or not _is_name(linear_kw.value, "izl") or quadratic_kw is None or not _is_name(quadratic_kw.value, "izq"):
        raise _Failure("C3 raw W return operand mapping changed")
    return {
        "result_assignment_ast_sha256": _ast_sha(result_nodes[0]),
        "generic_return_assignment_ast_sha256": _ast_sha(generic_kw),
        "sectors_return_assignment_ast_sha256": _ast_sha(sectors_kw),
        "inertial_axial_key": "actual_generic_source_numerators.inertial_axial",
        "sector_keys": ["full_signed_inertial_sectors_axial4.axial_linear", "full_signed_inertial_sectors_axial4.axial_quadratic"],
        "raw_W_not_recovered_as_E_times_P2": True,
    }


def _c3_formula(sources: Mapping[str, Any]) -> dict[str, Any]:
    recovery = sources[C3_RECOVERY]
    reference = sources[C3_REFERENCE]
    long = sources[C3_LONG]
    recovery_fn = _find(recovery["tree"], "recover_inputs")
    reference_fn = _find(reference["tree"], "quotient_cell")
    long_fn = _find(long["tree"], "source_quotients")
    if recovery_fn is None or reference_fn is None or long_fn is None:
        raise _Failure("C3 formula creators missing")
    history_assignment = _exact_histories_assignment(recovery_fn)
    pstar_nodes = [item for item in _assignment_nodes(recovery_fn, "S") if isinstance(item.value, ast.Attribute) and _attr_chain(item.value) == ["op", "Pstar"]]
    delta_nodes = [item for item in _assignment_nodes(recovery_fn, "de") if isinstance(item.value, ast.Attribute) and _attr_chain(item.value) == ["op", "reference", "delta"]]
    radius_nodes = _assignment_nodes(recovery_fn, "R")
    iz_nodes = _assignment_nodes(recovery_fn, "Iz")
    if len(pstar_nodes) != 1 or len(delta_nodes) != 1 or len(radius_nodes) != 1 or len(iz_nodes) != 1:
        raise _Failure("C3 Pstar/delta/radius/Iz assignments changed")
    raw_return = _c3_provider_return(recovery_fn)
    before_nodes = _assignment_nodes(long_fn, "before")
    p2_nodes = _assignment_nodes(long_fn, "p2")
    if len(before_nodes) != 1 or len(p2_nodes) != 1:
        raise _Failure("C3 source_quotients assignments changed")
    before = before_nodes[0].value
    if not (_is_call_chain(before, ["parameters", "add"]) and len(before.args) == 3 and _is_name(before.args[0], "f") and isinstance(before.args[2], ast.Call) and _is_call_chain(before.args[2], ["parameters", "scale"]) and len(before.args[2].args) == 2 and isinstance(before.args[2].args[1], ast.Attribute) and _attr_chain(before.args[2].args[1]) == ["op", "Pstar"]):
        raise _Failure("C3 Pstar inertial source assignment changed")
    quotient_calls = [
        item
        for item in ast.walk(long_fn)
        if isinstance(item, ast.Call)
        and _attr_chain(item.func) == ["parameters", "quotient"]
        and len(item.args) >= 3
        and _is_name(item.args[0], "f")
        and _is_name(item.args[2], "E")
    ]
    if len(quotient_calls) != 1:
        raise _Failure("C3 source_quotients E quotient route changed")
    p2_value = p2_nodes[0].value
    if not (_is_call_chain(p2_value, ["parameters", "scale"]) and len(p2_value.args) == 2 and _subscript(p2_value.args[0], "quotients", "p2") and isinstance(p2_value.args[1], ast.Attribute) and _attr_chain(p2_value.args[1]) == ["op", "source_radius"]):
        raise _Failure("C3 physical-radius P2 route changed")
    reference_p2 = _assignment_nodes(reference_fn, "p2")
    if len(reference_p2) != 1 or not (
        _is_call_chain(reference_p2[0].value, ["parameters", "scale"])
        and len(reference_p2[0].value.args) == 2
        and _subscript(reference_p2[0].value.args[0], "inertial", "p2")
        and isinstance(reference_p2[0].value.args[1], ast.Attribute)
        and _attr_chain(reference_p2[0].value.args[1]) == ["proxy", "source_radius"]
    ):
        raise _Failure("C3 reference physical-radius P2 route changed")
    return {
        "raw_W_source": {
            "source_path": recovery["path"],
            "source_sha": recovery["sha"],
            "function": "recover_inputs",
            "function_ast_sha256": _ast_sha(recovery_fn),
            "Iz_assignment_ast_sha256": _ast_sha(iz_nodes[0]),
            "Iz_assignment_source": _text(iz_nodes[0]),
            "provider_return": raw_return,
            "radius_assignment_ast_sha256": _ast_sha(radius_nodes[0]),
            "Pstar_assignment_ast_sha256": _ast_sha(pstar_nodes[0]),
            "delta_assignment_ast_sha256": _ast_sha(delta_nodes[0]),
            "history_assignment_ast_sha256": _ast_sha(history_assignment),
            "invS_applied_once": True,
        },
        "source_quotients_route": {
            "source_path": long["path"],
            "source_sha": long["sha"],
            "function": "source_quotients",
            "function_ast_sha256": _ast_sha(long_fn),
            "Pstar_assignment_ast_sha256": _ast_sha(before_nodes[0]),
            "E_quotient_call_ast_sha256": _ast_sha(quotient_calls[0]),
            "P2_assignment_ast_sha256": _ast_sha(p2_nodes[0]),
            "P2_assignment_is_not_raw_W": True,
        },
        "reference_provider": {
            "source_path": reference["path"],
            "source_sha": reference["sha"],
            "function": "quotient_cell",
            "function_ast_sha256": _ast_sha(reference_fn),
            "P2_assignment_ast_sha256": _ast_sha(reference_p2[0]),
            "P2_source_route_present": True,
        },
        "W_side": "R*(I_z_linear + Pstar*I_z_quadratic)",
        "P2_side": "W_side/E",
    }


def _route_call(fn: ast.AST, chain: list[str], args: list[tuple[str, Any]]) -> ast.Call | None:
    for item in ast.walk(fn):
        if isinstance(item, ast.Call) and _attr_chain(item.func) == chain and len(item.args) >= len(args):
            okay = True
            for actual, (kind, expected) in zip(item.args, args):
                if kind == "name" and not _is_name(actual, expected):
                    okay = False
                if kind == "subscript" and not _subscript(actual, expected[0], expected[1]):
                    okay = False
            if okay:
                return item
    return None


def _c3_routes(sources: Mapping[str, Any]) -> dict[str, Any]:
    original = _find(sources[C3_O3]["tree"], "query", "OriginalO3RcFiniteN")
    whole = _find(sources[C3_WHOLE_O3]["tree"], "query", "WholeZO3RcFiniteN")
    leading = _find(sources[C3_WHOLE]["tree"], "leading_packet", "WholeZAllNOuterRcFunctions")
    repair_leading = _find(sources[C3_REPAIR]["tree"], "leading_C3")
    repair_build = _find(sources[C3_REPAIR]["tree"], "build")
    compile_recovery = _find(sources[C3_LONG]["tree"], "compile_recovery")
    if any(item is None for item in (original, whole, leading, repair_leading, repair_build, compile_recovery)):
        raise _Failure("C3 route owner/source function missing")
    route_views = []
    for owner, fn, path in (
        ("OriginalO3RcFiniteN.query", original, C3_O3),
        ("WholeZO3RcFiniteN.query", whole, C3_WHOLE_O3),
    ):
        recover = _route_call(fn, ["reference", "long", "RECOVER"], [("name", "proxy"), ("subscript", ("source", "raw"))])
        quotient = _route_call(fn, ["original", "O3_QUOTIENTS"], []) or _route_call(fn, ["O3_QUOTIENTS"], [])
        if recover is None or quotient is None:
            raise _Failure(f"C3 O3 route shape changed: {owner}")
        route_views.append(
            {
                "owner": owner,
                "source_path": path,
                "source_sha": sources[path]["sha"],
                "function_ast_sha256": _ast_sha(fn),
                "recovery_call_ast_sha256": _ast_sha(recover),
                "quotient_call_ast_sha256": _ast_sha(quotient),
            }
        )
    frontend_calls = [
        item
        for item in ast.walk(leading)
        if isinstance(item, ast.Call)
        and isinstance(item.func, ast.Subscript)
        and _is_name(item.func.value, "FRONTENDS")
    ]
    if len(frontend_calls) != 2:
        raise _Failure("whole-Z leading FRONTENDS callback shape changed")
    p0_nodes = [item for item in _assignment_nodes(repair_leading, "P0") if _is_call_chain(item.value, ["third_leaf"]) and len(item.value.args) == 1 and _subscript(item.value.args[0], "prior", "independent_P0")]
    histories_nodes = _assignment_nodes(repair_leading, "histories")
    complete_nodes = _assignment_nodes(repair_build, "complete")
    if len(p0_nodes) != 1 or len(histories_nodes) != 1 or len(complete_nodes) != 1:
        raise _Failure("C3 repair-band P0/history producer shape changed")
    return {
        "O3_route_callbacks": route_views,
        "compiled_recovery_route": {
            "source_path": C3_LONG,
            "source_sha": sources[C3_LONG]["sha"],
            "function": "compile_recovery",
            "function_ast_sha256": _ast_sha(compile_recovery),
        },
        "whole_Z_leading_packet": {
            "owner": "WholeZAllNOuterRcFunctions.leading_packet",
            "source_path": C3_WHOLE,
            "source_sha": sources[C3_WHOLE]["sha"],
            "function_ast_sha256": _ast_sha(leading),
            "callback_ast_sha256": [_ast_sha(item) for item in frontend_calls],
        },
        "C3_repair_band_producer": {
            "source_path": C3_REPAIR,
            "source_sha": sources[C3_REPAIR]["sha"],
            "leading_function": "leading_C3",
            "leading_function_ast_sha256": _ast_sha(repair_leading),
            "P0_assignment_ast_sha256": _ast_sha(p0_nodes[0]),
            "history_assignment_ast_sha256": _ast_sha(histories_nodes[0]),
            "build_function": "build",
            "build_function_ast_sha256": _ast_sha(repair_build),
            "complete_history_assignment_ast_sha256": _ast_sha(complete_nodes[0]),
        },
    }


def _pressure(sources: Mapping[str, Any]) -> dict[str, Any]:
    background = sources[NATIVE_BACKGROUND]
    recovery = sources[C3_RECOVERY]
    init = _find(background["tree"], "__init__", "CurrentOriginalTrueOuterC6C5Functions")
    transport = _find(background["tree"], "transport_geometry")
    recovery_fn = _find(recovery["tree"], "recover_inputs")
    if init is None or transport is None or recovery_fn is None:
        raise _Failure("pressure route creators missing")
    native_nodes = [item for item in _subscript_assignment_nodes(init, "fields") if isinstance(item.value, ast.Call) and _is_call_chain(item.value, ["a", "add"]) and len(item.value.args) == 2 and _subscript(item.targets[0], "fields", "pressure") and _subscript(item.value.args[0], "fields", "P0") and _subscript(item.value.args[1], "fields", "p")]
    c3_nodes = [item for item in _assignment_nodes(recovery_fn, "pressure") if isinstance(item.value, ast.Call) and _is_name(item.value.func, "add") and len(item.value.args) == 2 and isinstance(item.value.args[0], ast.Attribute) and _attr_chain(item.value.args[0]) == ["op", "P0"] and _is_name(item.value.args[1], "p")]
    if len(native_nodes) != 1 or len(c3_nodes) != 1 or not {"incoming", "local", "memory"}.issubset(_names(transport)):
        raise _Failure("pressure restoration AST route changed")
    return {
        "native": {"source_path": background["path"], "source_sha": background["sha"], "function": "CurrentOriginalTrueOuterC6C5Functions.__init__", "assignment_ast_sha256": _ast_sha(native_nodes[0]), "transport_function_ast_sha256": _ast_sha(transport)},
        "rebuilt_C3": {"source_path": recovery["path"], "source_sha": recovery["sha"], "function": "recover_inputs", "assignment_ast_sha256": _ast_sha(c3_nodes[0]), "independent_P0": True},
        "pressure_side": "P0 + p",
    }


def _projection(sources: Mapping[str, Any], c3_literals: Mapping[str, Any]) -> dict[str, Any]:
    c2 = _find(sources[C2_PROJECTION]["tree"], "source_packet", "CurrentC2PhasePrimitives")
    c3 = _find(sources[C3_PROJECTION]["tree"], "source_packet", "CurrentC3SourceTargets")
    if c2 is None or c3 is None:
        raise _Failure("projection creators missing")
    def one(fn: ast.AST, factors: dict[int, int]) -> ast.Assign:
        rows = _subscript_assignment_nodes(fn, "projected")
        if len(rows) != 1 or not isinstance(rows[0].value, ast.List):
            raise _Failure("projection assignment multiplicity changed")
        items = rows[0].value.elts
        for index, factor in factors.items():
            if index >= len(items) or not _factor_row(items[index], index, factor):
                raise _Failure("projection factor AST changed")
        return rows[0]
    c2_assignment = one(c2, {2: 2})
    c3_assignment = one(c3, {2: 2, 3: 6})
    if sum(_factor_row(item, 2, 2) for item in ast.walk(c2)) != 1 or sum(_factor_row(item, 2, 2) for item in ast.walk(c3)) != 1 or sum(_factor_row(item, 3, 6) for item in ast.walk(c3)) != 1:
        raise _Failure("ordinary-Z factorial application is not exactly once")
    return {
        "C2": {"source_path": C2_PROJECTION, "source_sha": sources[C2_PROJECTION]["sha"], "function": "CurrentC2PhasePrimitives.source_packet", "function_ast_sha256": _ast_sha(c2), "assignment_ast_sha256": _ast_sha(c2_assignment), "factorial_applied_once": {"2": True}},
        "C3": {"source_path": C3_PROJECTION, "source_sha": sources[C3_PROJECTION]["sha"], "function": "CurrentC3SourceTargets.source_packet", "function_ast_sha256": _ast_sha(c3), "assignment_ast_sha256": _ast_sha(c3_assignment), "factorial_applied_once": {"2": True, "3": True}},
        "factorial_map": {str(key): value for key, value in FACTORIALS.items()},
        "accepted_recipe_rows": c3_literals["p2_root_selectors"],
    }


def _factor_row(node: Any, row_index: int, factor: int) -> bool:
    if not isinstance(node, ast.BinOp) or not isinstance(node.op, ast.Mult):
        return False
    pairs = ((node.left, node.right), (node.right, node.left))
    return any(_subscript(first, "row", row_index) and _is_const(second, factor) for first, second in pairs)


def _views(native: Mapping[str, Any], c3: Mapping[str, Any], pressure: Mapping[str, Any]) -> dict[str, Any]:
    def side(name: str, formula: Mapping[str, Any], literals: Mapping[str, Any], raw_w: Mapping[str, Any]) -> dict[str, Any]:
        rows = literals["p2_root_selectors"]
        raw_views = []
        for row in rows:
            raw_views.append(
                {
                    "order": row["order"],
                    "derivative_order": row["order"],
                    "dependency_order": row["dependency_order"],
                    "raw_W_source": raw_w,
                    "source_path": raw_w.get("source_path"),
                    "source_sha": raw_w.get("source_sha"),
                    "source_function_ast_sha256": raw_w.get("function_ast_sha256"),
                    "source_assignment_ast_sha256": raw_w.get(
                        "assignment_ast_sha256", raw_w.get("Iz_assignment_ast_sha256")
                    ),
                    "p2_root_selector": row,
                    "p2_numerator_is_raw_W": False,
                    "derivative_creator_is_actual_source": True,
                    "view_kind": "typed source-assignment derivative of complete W_side",
                }
            )
        return {
            "side": name,
            "W_side": "R*(I_z_linear + Pstar*I_z_quadratic)",
            "P2_side": "W_side/E",
            "side_local_factors": ["R", "Pstar", "E", "delta", "P0", "p"],
            "raw_W_assignment_source": raw_w,
            "raw_W_derivative_views": raw_views,
            "p2_quotient_projection": {"rows": rows, "is_raw_W": False},
            "pressure_route": pressure[name],
        }
    return {"native": side("native", native["formula"], native["literal_handles"], native["raw_W_source"]), "rebuilt_C3": side("rebuilt_C3", c3["formula"], c3["literal_handles"], c3["raw_W_source"])}


def _audit_zero() -> dict[str, Any]:
    return {
        "new_mathematical_proofs": 0,
        "new_algebra_family_calls": 0,
        "numeric_calls": 0,
        "source_provider_calls": 0,
        "old_math_calls": 0,
        "checker_calls": 0,
        "normalizer_calls": 0,
        "reducer_calls": 0,
        "quadrature_calls": 0,
        "recurrence_calls": 0,
        "field_evaluations": 0,
        "integral_evaluations": 0,
        "provider_calls": 0,
    }


def _counts(accepted: bool) -> dict[str, int]:
    result = {key: 0 for key in COUNT_KEYS}
    if accepted:
        result.update(
            source_port_rows=8,
            native_source_port_rows=4,
            C3_source_port_rows=4,
            projection_rows=8,
            pressure_route_sides=2,
        )
    return result


def _failure(events: list[dict[str, Any]], error: str, outer: Mapping[str, Any] | None) -> dict[str, Any]:
    flags = {key: False for key in NEW_FLAGS}
    result = {
        "schema": RESULT_SCHEMA,
        "component": COMPONENT,
        "status": "failed",
        "passed": False,
        "source_only": True,
        "metadata_only": False,
        "runtime_model_metadata": None,
        "events": list(events),
        "journal_events_prefix": list(events),
        "failure_prefix_events": list(events),
        "blockers": [error],
        "new_flags": flags,
        "flags": dict(flags),
        **flags,
        "actual_counts": _counts(False),
        "audit_counts": _counts(False),
        "source_port_rows": 0,
        "native_source_port_rows": 0,
        "C3_source_port_rows": 0,
        "projection_rows": 0,
        "pressure_route_sides": 0,
        "actual_projection_rows": 0,
        "actual_pressure_route_sides": 0,
        "new_mathematical_proofs": 0,
        "new_algebra_family_calls": 0,
        "numeric_calls": 0,
        "source_provider_calls": 0,
        "old_math_calls": 0,
        "actual_projection_orders": [],
        "audit": _audit_zero(),
        "source_packet_defined": False,
        "source_packet_definition_accepted": False,
        "source_packet_accepted_owned": False,
        "source_packet_accepted_inherited": False,
        "full_goal_complete": False,
        "numeric_values_computed": False,
        "numeric_field_values_present": False,
    }
    if isinstance(outer, Mapping):
        result.update(input_digest=outer.get("input_digest"), input_raw_sha256=outer.get("input_raw_sha256"), parent_result_raw_sha256=outer.get("parent_result_raw_sha256"))
    return result


def derive(source: Mapping[str, Any], journal: Callable[[dict[str, Any]], Any] | None = None) -> dict[str, Any]:
    """Bind one frozen outer wrapper; reject raw v1 packets passed directly."""

    events: list[dict[str, Any]] = []
    outer = source if isinstance(source, Mapping) else None
    def emit(name: str, **payload: Any) -> None:
        _emit(events, journal, {"event": name, **payload})
    try:
        if not isinstance(source, Mapping) or source.get("schema") != INPUT_SCHEMA:
            raise _Failure("derive requires frozen outer source wrapper")
        payload = source.get("input")
        if not isinstance(payload, Mapping) or payload is source or payload.get("schema") != INPUT_SCHEMA:
            raise _Failure("derive requires exactly one source['input'] v1 payload")
        trace = _mapping(payload.get("trace"), "trace")
        if trace.get("sha") != TRACE_SHA:
            raise _Failure("trace SHA pin changed")
        emit("outer_wrapper_unwrapped_once", input_schema=INPUT_SCHEMA)
        records = _mapping(payload.get("creator_source_records"), "creator_source_records")
        sources: dict[str, Any] = {}
        for path in REQUIRED_SOURCES:
            if path not in records:
                raise _Failure(f"required source record missing: {path}")
            sources[path] = _bind_source(records[path], path)
            emit("creator_source_bound", path=path, source_sha256=sources[path]["sha"], module_ast_sha256=sources[path]["module_ast_sha256"])
        slices = _mapping(payload.get("literal_slices"), "literal_slices")
        native_literals = _native_literals(_mapping(slices.get("native"), "literal_slices.native"))
        emit("native_P2_quotient_selectors_bound", root_ids=list(NATIVE_ROOT_IDS), coordinate_node=241536)
        c3_literals = _c3_literals(_mapping(slices.get("C3"), "literal_slices.C3"))
        emit("rebuilt_C3_P2_quotient_selectors_bound", recipe_ids=list(C3_ROOT_IDS), coordinate_node=4582)
        native_formula = _native_formula(sources)
        c3_formula = _c3_formula(sources)
        routes = _c3_routes(sources)
        pressure = _pressure(sources)
        projection = _projection(sources, c3_literals)
        emit("whole_function_raw_W_sources_bound", native_assignment=native_formula["assignment_ast_sha256"], c3_assignment=c3_formula["raw_W_source"]["Iz_assignment_ast_sha256"])
        emit("pressure_restoration_route_bound", sides=2, pressure_side="P0 + p")
        emit("ordinary_Z_projection_bound", orders=list(ORDINARY_Z), factorials={str(k): v for k, v in FACTORIALS.items()})
        native_ref = {"formula": native_formula, "literal_handles": native_literals, "raw_W_source": native_formula, "ordinary_Z_projection": native_literals["p2_root_selectors"]}
        c3_ref = {"formula": c3_formula, "literal_handles": c3_literals, "raw_W_source": c3_formula["raw_W_source"], "ordinary_Z_projection": c3_literals["p2_root_selectors"]}
        flags = {key: True for key in NEW_FLAGS}
        counts = _counts(True)
        audit = _audit_zero()
        views = _views(native_ref, c3_ref, pressure)
        emit("source_binding_accepted", flags=flags, actual_counts=counts, projection_rows=8)
        return {
            "schema": RESULT_SCHEMA,
            "component": COMPONENT,
            "status": "passed",
            "passed": True,
            "source_only": True,
            "metadata_only": False,
            "runtime_model_metadata": None,
            "events": events,
            "journal_events": list(events),
            "new_flags": flags,
            "flags": dict(flags),
            **flags,
            "actual_counts": counts,
            "audit_counts": dict(counts),
            "source_port_rows": counts["source_port_rows"],
            "native_source_port_rows": counts["native_source_port_rows"],
            "C3_source_port_rows": counts["C3_source_port_rows"],
            "projection_rows": counts["projection_rows"],
            "pressure_route_sides": counts["pressure_route_sides"],
            "actual_projection_rows": counts["projection_rows"],
            "actual_pressure_route_sides": counts["pressure_route_sides"],
            "new_mathematical_proofs": 0,
            "new_algebra_family_calls": 0,
            "numeric_calls": 0,
            "source_provider_calls": 0,
            "old_math_calls": 0,
            "actual_projection_orders": list(ORDINARY_Z),
            "audit": audit,
            "typed_source_references": {"native": native_ref, "rebuilt_C3": c3_ref},
            "source_selector_ast_audit": {path: sources[path]["selector_ast_audit"] for path in REQUIRED_SOURCES},
            "whole_function_pressure_numerator_views": views,
            "C3_source_routes": routes,
            "projection": projection,
            "dependency_order": "j+1 for each ordinary Z order j",
            "C3_history_normalization": "invS m/k conditional assignment bound once",
            "source_packet_defined": True,
            "source_packet_definition_accepted": True,
            "source_packet_accepted_owned": True,
            "source_packet_accepted_inherited": False,
            "full_goal_complete": False,
            "numeric_values_computed": False,
            "numeric_field_values_present": False,
            "input_digest": source.get("input_digest"),
            "input_raw_sha256": source.get("input_raw_sha256"),
            "parent_result_raw_sha256": source.get("parent_result_raw_sha256"),
            "unresolved": [
                "No native/C3 equality of E, phase, radius, Pstar, delta, P0, or p is asserted.",
                "No P2 value, pressure closure, bound, or signed gap zero is asserted.",
            ],
        }
    except _Failure as exc:
        return _failure(events, str(exc), outer)
    except Exception as exc:
        return _failure(events, f"unexpected source-binding failure: {type(exc).__name__}: {exc}", outer)

